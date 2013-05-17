//============================================================================
//                                  I B E X
// Linear Relaxation Contractor
// File        : ibex_CtcLinearRelaxation.cpp
// Author      : Bertrand Neveu , Gilles Trombettoni, Jordan Ninin
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Nov 14, 2012
// Last Update : May 15, 2013
//============================================================================

#include "ibex_CtcLinearRelaxation.h"

namespace ibex {


const double default_max_diam_box=1e6;
const int default_max_time_out=100;
const int default_max_iter=100;


CtcLinearRelaxation::CtcLinearRelaxation(const System& sys, int goal_ctr,
		ctc_mode cmode, int max_iter, int time_out, double max_diam_box):
				  Ctc(sys.nb_var), sys(sys), goal_ctr(goal_ctr),
				  cmode(cmode),  max_iter(max_iter),
				  max_time_out(time_out), max_diam_box(max_diam_box)
{

	mylinearsolver = new LinearSolver(sys.nb_var, sys.nb_ctr, max_iter, time_out);

	// ============================================================

	mylinearsolver.setSense(LinearSolver::MINIMIZE);
	mylinearsolver.setMaxIter(max_iter);
	mylinearsolver.setMaxTimeOut(time_out);
	mylinearsolver.setEpsilon(1e-10);

	primal_solution = new double[sys.nb_var];

}

CtcLinearRelaxation::~CtcLinearRelaxation () {
	delete mylinearsolver;
	delete[] primal_solution;
}


void CtcLinearRelaxation::optimizer(IntervalVector & box){

	Interval opt(0);
	int* inf_bound = new int[sys.nb_var]; // indicator inf_bound = 1 means the inf bound is feasible or already contracted, call to simplex useless (cf Baharev)
	int* sup_bound = new int[sys.nb_var]; // indicator sup_bound = 1 means the sup bound is feasible or already contracted, call to simplex useless

	for (int i=0; i<sys.nb_var ;i++) {inf_bound[i]=0;sup_bound[i]=0;}
	if (!sys.goal)   sup_bound[sys.nb_var-1]=1;  // in case of optimization, for y the left bound only is contracted

	int firsti=0;
	// in the case of lower_bounding, only the left bound of y is contracted
	if ((!sys.goal) && cmode==ONLY_Y) firsti=2*sys.nb_var-1;

	int nexti=-1;  // the next variable to be contracted
	int infnexti=0; // the bound to be contracted contract  infnexti=0 for the lower bound, infnexti=1 for the upper bound


	for(int ii=firsti;ii<(2*sys.nb_var);ii++){  // at most 2*n calls
		int i= ii/2;
		if (nexti != -1) i=nexti;
		//    cout << " i "<< i << " infnexti " << infnexti << " infbound " << inf_bound[i] << " supbound " << sup_bound[i] << endl;

		if (infnexti==0 && inf_bound[i]==0)  // computing the left bound : minimizing x_i
		{
			inf_bound[i]=1;
			LinearSolver::Status_Sol stat = run_simplex(box, LinearSolver::MINIMIZE, i, sys.nb_var, opt,box[i].lb());

			if( stat == LinearSolver::OPTIMAL )
			{
				if(opt.lb()>box[i].ub())  throw EmptyBoxException();

				choose_next_variable(box,nexti,infnexti, inf_bound, sup_bound);

				if(opt.lb() > box[i].lb() ){
					box[i]=Interval(opt.lb(),box[i].ub());
					mylinearsolver.setBoundVar(i,box[i]);
				}
			}
			else if (stat == LinearSolver::INFEASIBLE)
			{
				// the infeasibility is proved, the EmptyBox exception is raised
				throw EmptyBoxException();
			}
			else if (stat == LinearSolver::UNKNOWN)
			{
				int next=-1;
				for (int j=0;j<sys.nb_var;j++)
				{
					if (inf_bound[j]==0) {
						nexti=j;  next=0;  infnexti=0;
						break;
					}
					else if  (sup_bound[j]==0) {
						nexti=j;  next=0;  infnexti=1;
						break;
					}
				}
				if (next==-1) {
					break;
				}
			}

		}
		else if(infnexti==1 && sup_bound[i]==0)// computing the right bound :  maximizing x_i
		{
			sup_bound[i]=1;
			LinearSolver::Status_Sol stat= run_simplex(box, LinearSolver::MAXIMIZE, i, sys.nb_var, opt, box[i].ub());

			if( stat == LinearSolver::OPTIMAL ){
				if(opt.ub() <box[i].lb())   throw EmptyBoxException();

				choose_next_variable(box,nexti,infnexti, inf_bound, sup_bound);

				if (opt.ub() < box[i].ub()) {
					box[i] =Interval( box[i].lb(), opt.ub());
					mylinearsolver.setBoundVar(i,box[i]);
				}
			}
			else if(stat == LinearSolver::INFEASIBLE)
			{
				// the infeasibility is proved,  the EmptyBox exception is raised
				throw EmptyBoxException();
			}

			else if (stat == LinearSolver::UNKNOWN)
			{
				int next=-1;
				for (int j=0;j<sys.nb_var;j++)
				{if (inf_bound[j]==0) {
					nexti=j; next=0;infnexti=0;
					break;
				}
				else if  (sup_bound[j]==0) {
					nexti=j; next=0;infnexti=1;
					break;
				}
				}
				if (next==-1) break;
			}
		}
		else {
			break;
		}  // no more call to the LP Solver
	}

	delete [] inf_bound;
	delete [] sup_bound;

}





LinearSolver::Status_Sol CtcLinearRelaxation::run_simplex(IntervalVector& box, LinearSolver::Sense sense, int var, int n,
		Interval& obj, double bound){

	//  mylinearsolver.writeFile("dump.lp");
	int nr=mylinearsolver.getNbRows();

	if(sense==LinearSolver::MINIMIZE)
		mylinearsolver.setVarObj(var, 1.0);
	else
		mylinearsolver.setVarObj(var, -1.0);

	LinearSolver::Status_Sol stat;

	//    mylinearsolver.writeFile("dump.lp");
	//    system("cat dump.lp");
	stat = mylinearsolver.solve();


	if(stat == LinearSolver::OPTIMAL){
		if( (sense==LinearSolver::MINIMIZE && mylinearsolver.getObjValue()<=bound) ||
			(sense==LinearSolver::MAXIMIZE && -mylinearsolver.getObjValue()>=bound) )
		{
			stat = LinearSolver::UNKNOWN;
		}
	}

	// Neumaier - Shcherbina postprocessing
	if(stat == LinearSolver::OPTIMAL){

		// the dual solution : used to compute the bound
		double * dual_solution = new double[mylinearsolver.getNbRows()];
		LinearSolver::Status stat_dual = mylinearsolver.getDualSol(dual_solution);

		Matrix As (mylinearsolver.getNbRows(),sys.nb_var) ;
		LinearSolver::Status stat_As = mylinearsolver.getCoefConstraint(As);

		IntervalVector B(mylinearsolver.getNbRows());
		LinearSolver::Status stat_B = mylinearsolver.getB(B);

		if (stat_dual==LinearSolver::OK && stat_As==LinearSolver::OK && stat_B==LinearSolver::OK)
			NeumaierShcherbina_postprocessing(n, nr, var, obj, box, As, B, (sense==LinearSolver::MINIMIZE),dual_solution);
		else
			stat = LinearSolver::INFEASIBLE;

		delete [] dual_solution;
	}

	// infeasibility test  cf Neumaier Shcherbina paper
	if(stat == LinearSolver::INFEASIBLE)
	{
		stat = LinearSolver::UNKNOWN;

		double * infeasible_dir= new double[mylinearsolver.getNbRows()];
		LinearSolver::Status stat1 = mylinearsolver.getInfeasibleDir(infeasible_dir);

		Matrix As (mylinearsolver.getNbRows(),sys.nb_var) ;
		LinearSolver::Status stat2 = mylinearsolver.getCoefConstraint(As);

		IntervalVector B(mylinearsolver.getNbRows());
		LinearSolver::Status stat3 = mylinearsolver.getB(B);

		if (stat1==LinearSolver::OK && stat2==LinearSolver::OK && stat3==LinearSolver::OK &&
				NeumaierShcherbina_infeasibilitytest (n, nr, box, As, B, infeasible_dir))
		{
			stat = LinearSolver::INFEASIBLE;
		}

		delete[] infeasible_dir;
	}


	mylinearsolver.setVarObj(var, 0.0);

	return stat;

}




// Neumaier Shcherbina postprocessing in case of optimal solution found : the result obj is made reliable
void CtcLinearRelaxation::NeumaierShcherbina_postprocessing (int n, int nr, int var, Interval & obj, IntervalVector& box,
		Matrix & As, IntervalVector& B, bool minimization, double * dual_solution) {

	IntervalVector Rest(n);
	for (int i =0; i< n ; i++)
		Rest[i] = As.row(i) * dual_solution ; // Rest = Transpose(As)*Lambda

	Rest[var] += ( (minimization)? -1 : 1); // because C is a vector of zero except for the coef "var"

	if(minimization)
		obj = dual_solution * B - Rest * box;
	else
		obj = -(dual_solution * B - Rest * box);
}


// Neumaier Shcherbina postprocessing in case of infeasibilty found by LP  returns true if the infeasibility is proved
bool  CtcLinearRelaxation::NeumaierShcherbina_infeasibilitytest (int n, int nr, IntervalVector& box, Matrix & As, IntervalVector& B, double * infeasible_dir)
{

	IntervalVector Rest(n);
	for (int i =0; i< n ; i++)
		Rest[i] = As.row(i) * infeasible_dir ;  // Rest = Transpose(As) * Lambda

	Interval d= Rest *box - infeasible_dir * B;

	// if 0 does not belong to d, the infeasibility is proved
	if (d.lb() > 0 || d.ub() <0)
		return true;
	else
		return false;
}


// check if the constraint is satisfied in the box : in this case, no linear relaxation is made.
bool CtcLinearRelaxation::isInner(IntervalVector & box,const System& sys, int j){
	Interval eval=sys.ctrs[j].f.eval(box);

	if((sys.ctrs[j].op==LEQ && eval.ub() > 0) || (sys.ctrs[j].op==LT && eval.ub() >= 0) ||
			(sys.ctrs[j].op==GEQ && eval.lb() < 0) || (sys.ctrs[j].op==GT && eval.lb() <= 0))
		return false;
	else
		return true;
}



// The Achterberg heuristic for choosing the next variable (nexti) and its bound (infnexti) to be contracted (cf Baharev paper)
// and updating the indicators if a bound has been found feasible (with the precision prec_bound)
// called only when a solution is found by the LP solver (use of primal_solution)

void CtcLinearRelaxation::choose_next_variable ( IntervalVector & box, int & nexti, int & infnexti, int* inf_bound, int* sup_bound)
{

	// the primal solution : used by choose_next_variable
	LinearSolver::Status stat_prim = mylinearsolver.getPrimalSol(primal_solution);
	
	if (stat_prim==OK) {
		double prec_bound = 1.e-8; // relative precision for the indicators
		double delta=1.e100;
		double deltaj=delta;

		for (int j=0;j<sys.nb_var;j++)	{

			if (inf_bound[j]==0) {
				deltaj= fabs (primal_solution[j]- box[j].lb());
				if ((fabs (box[j].lb()) < 1 && deltaj < prec_bound) ||
					(fabs (box[j].lb()) >= 1 && fabs (deltaj /(box[j].lb())) < prec_bound))	{
					inf_bound[j]=1;
				}
				if (inf_bound[j]==0 && deltaj < delta) 	{
					nexti=j; infnexti=0;delta=deltaj;
				}
			}

			if (sup_bound[j]==0) {
				deltaj = fabs (primal_solution[j]- box[j].ub());

				if ((fabs (box[j].ub()) < 1 && deltaj < prec_bound) 	||
					(fabs (box[j].ub()) >= 1 && fabs (deltaj/(box[j].ub())) < prec_bound)) {
					sup_bound[j]=1;
				}
				if (sup_bound[j]==0 && deltaj < delta) {
					nexti=j; infnexti=1;delta=deltaj;
				}
			}

		}
	}
	else {
		// TODO find something to do by default 
	}

}


}// end namespace ibex

