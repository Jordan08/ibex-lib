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


CtcLinearRelaxation::CtcLinearRelaxation(const System& sys1,
		ctc_mode cmode, int max_iter, int time_out, double eps, double max_diam_box1):
		 Ctc(sys1.nb_var),
		 sys(sys1),
		 cmode(cmode),
		 max_diam_box(max_diam_box1)
//		 , primal_solution(sys.nb_var)
//		 , mylinearsolver(sys.nb_var, sys.nb_ctr, max_iter, time_out, eps)
{

	goal_ctr= (sys1.goal!=NULL) ? 0 : -1 ;  // TODO to modify with the value store in System
	goal_var= (sys1.goal!=NULL) ? (sys1.nb_var-1) : -1 ;  // TODO to modify with the value store in System
	mylinearsolver = new LinearSolver(sys1.nb_var, sys1.nb_ctr, max_iter, time_out, eps);
}

CtcLinearRelaxation::~CtcLinearRelaxation () {
	delete mylinearsolver;
}


void CtcLinearRelaxation::optimizer(IntervalVector & box){

	Interval opt(0.0);
	int* inf_bound = new int[sys.nb_var]; // indicator inf_bound = 1 means the inf bound is feasible or already contracted, call to simplex useless (cf Baharev)
	int* sup_bound = new int[sys.nb_var]; // indicator sup_bound = 1 means the sup bound is feasible or already contracted, call to simplex useless

	for (int i=0; i<sys.nb_var ;i++) {inf_bound[i]=0;sup_bound[i]=0;}
	if (sys.goal!=NULL)   sup_bound[goal_var]=1;  // in case of optimization, for y the left bound only is contracted

	int firsti=0;
	// in the case of lower_bounding, only the left bound of y is contracted
	if ((sys.goal!=NULL) && cmode==ONLY_Y) firsti=2*sys.nb_var-1;

	int nexti=-1;  // the next variable to be contracted
	int infnexti=0; // the bound to be contracted contract  infnexti=0 for the lower bound, infnexti=1 for the upper bound


	for(int ii=firsti;ii<(2*sys.nb_var);ii++){  // at most 2*n calls
		int i= ii/2;
		if (nexti != -1) i=nexti;
//  cout << " i "<< i << " infnexti " << infnexti << " infbound " << inf_bound[i] << " supbound " << sup_bound[i] << endl;

		if (infnexti==0 && inf_bound[i]==0)  // computing the left bound : minimizing x_i
		{
			inf_bound[i]=1;
			LinearSolver::Status_Sol stat = run_simplex(box, LinearSolver::MINIMIZE, i, opt,box[i].lb());

			if( stat == LinearSolver::OPTIMAL )
			{
				if(opt.lb()>box[i].ub())  throw EmptyBoxException();

				if (!choose_next_variable(box,nexti,infnexti, inf_bound, sup_bound)) {
					break;
				}

				if(opt.lb() > box[i].lb() ){
					box[i]=Interval(opt.lb(),box[i].ub());
					mylinearsolver->setBoundVar(i,box[i]);
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
				if (next==-1)  break;
			}

		}
		else if (infnexti==1 && sup_bound[i]==0)// computing the right bound :  maximizing x_i
		{
			sup_bound[i]=1;
			LinearSolver::Status_Sol stat= run_simplex(box, LinearSolver::MAXIMIZE, i, opt, box[i].ub());

			if( stat == LinearSolver::OPTIMAL )
			{
				if(opt.ub() <box[i].lb())   throw EmptyBoxException();

				if (!choose_next_variable(box,nexti,infnexti, inf_bound, sup_bound)) {
					break;
				}

				if (opt.ub() < box[i].ub()) {
					box[i] =Interval( box[i].lb(), opt.ub());
					mylinearsolver->setBoundVar(i,box[i]);
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
				for (int j=0;j<sys.nb_var;j++)  {
					if (inf_bound[j]==0) {
						nexti=j;  next=0;  infnexti=0;
						break;
					}
					else if  (sup_bound[j]==0) {
						nexti=j;  next=0;  infnexti=1;
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



LinearSolver::Status_Sol CtcLinearRelaxation::run_simplex(IntervalVector& box,
		LinearSolver::Sense sense, int var, Interval& obj, double bound){


	if(sense==LinearSolver::MINIMIZE)
		mylinearsolver->setVarObj(var, 1.0);
	else
		mylinearsolver->setVarObj(var, -1.0);

	//    mylinearsolver->writeFile("dump.lp");
	//    system("cat dump.lp");
	LinearSolver::Status_Sol stat = mylinearsolver->solve();


	if(stat == LinearSolver::OPTIMAL){
		if( ((sense==LinearSolver::MINIMIZE) && (  mylinearsolver->getObjValue() <=bound)) ||
			((sense==LinearSolver::MAXIMIZE) && ((-mylinearsolver->getObjValue())>=bound)) )
		{
			stat = LinearSolver::UNKNOWN;
		}
	}

	// Neumaier - Shcherbina postprocessing
	if(stat == LinearSolver::OPTIMAL){

		// the dual solution : used to compute the bound
		Vector dual_solution(mylinearsolver->getNbRows());
		LinearSolver::Status stat_dual = mylinearsolver->getDualSol(dual_solution);

		Matrix A_trans (sys.nb_var,mylinearsolver->getNbRows()) ;
		LinearSolver::Status stat_A = mylinearsolver->getCoefConstraint_trans(A_trans);

		IntervalVector B(mylinearsolver->getNbRows());
		LinearSolver::Status stat_B = mylinearsolver->getB(B);

		if (stat_dual==LinearSolver::OK && stat_A==LinearSolver::OK && stat_B==LinearSolver::OK)
			NeumaierShcherbina_postprocessing( var, obj, box, A_trans, B, LinearSolver::MINIMIZE,dual_solution);
		else
			stat = LinearSolver::INFEASIBLE;


	}

	// infeasibility test  cf Neumaier Shcherbina paper
	if(stat == LinearSolver::INFEASIBLE)
	{
		stat = LinearSolver::UNKNOWN;

		Vector infeasible_dir(mylinearsolver->getNbRows());
		LinearSolver::Status stat1 = mylinearsolver->getInfeasibleDir(infeasible_dir);

		Matrix A_trans (sys.nb_var,mylinearsolver->getNbRows()) ;
		LinearSolver::Status stat2 = mylinearsolver->getCoefConstraint_trans(A_trans);

		IntervalVector B(mylinearsolver->getNbRows());
		LinearSolver::Status stat3 = mylinearsolver->getB(B);

		if (stat1==LinearSolver::OK && stat2==LinearSolver::OK && stat3==LinearSolver::OK &&
				NeumaierShcherbina_infeasibilitytest ( box, A_trans, B, infeasible_dir))
		{
			stat = LinearSolver::INFEASIBLE;
		}
	}


	mylinearsolver->setVarObj(var, 0.0);

	return stat;

}




// Neumaier Shcherbina postprocessing in case of optimal solution found : the result obj is made reliable
void CtcLinearRelaxation::NeumaierShcherbina_postprocessing ( int var, Interval & obj, IntervalVector& box,
		Matrix & A_trans, IntervalVector& B, LinearSolver::Sense minimization, Vector & dual_solution) {

	IntervalVector Rest(sys.nb_var);
	for (int i =0; i< sys.nb_var ; i++)
		Rest[i] = A_trans.row(i) * dual_solution ; // Rest = Transpose(A)*Lambda

	Rest[var] += ( (minimization==LinearSolver::MINIMIZE)? -1 : 1); // because C is a vector of zero except for the coef "var"

	if (minimization==LinearSolver::MINIMIZE)
		obj = dual_solution * B - Rest * box;
	else
		obj = -(dual_solution * B - Rest * box);
}


// Neumaier Shcherbina postprocessing in case of infeasibilty found by LP  returns true if the infeasibility is proved
bool  CtcLinearRelaxation::NeumaierShcherbina_infeasibilitytest ( IntervalVector& box,
			Matrix & A_trans, IntervalVector& B, Vector & infeasible_dir)
{

	IntervalVector Rest(sys.nb_var);
	for (int i =0; i< sys.nb_var; i++)
		Rest[i] = A_trans.row(i) * infeasible_dir ;  // Rest = Transpose(As) * Lambda

	Interval d= Rest *box - infeasible_dir * B;

	// if 0 does not belong to d, the infeasibility is proved
	if ((d.lb() > 0) ||( d.ub() <0))
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



bool CtcLinearRelaxation::choose_next_variable ( IntervalVector & box,
		int & nexti, int & infnexti, int* inf_bound, int* sup_bound)
{
	bool found = false;

	// the primal solution : used by choose_next_variable
	Vector primal_solution(sys.nb_var);
	LinearSolver::Status stat_prim = mylinearsolver->getPrimalSol(primal_solution);
	
	if (stat_prim==LinearSolver::OK) {
		// The Achterberg heuristic for choosing the next variable (nexti) and its bound (infnexti) to be contracted (cf Baharev paper)
		// and updating the indicators if a bound has been found feasible (with the precision prec_bound)
		// called only when a primal solution is found by the LP solver (use of primal_solution)

		double prec_bound = mylinearsolver->getEpsilon(); // relative precision for the indicators TODO change with the precision of the optimizer
		double delta=1.e100;
		double deltaj=delta;

		for (int j=0;j<sys.nb_var;j++)	{

			if (inf_bound[j]==0) {
				deltaj= fabs(primal_solution[j]- box[j].lb());
				if ((fabs (box[j].lb()) < 1 && deltaj < prec_bound) ||
					(fabs (box[j].lb()) >= 1 && fabs (deltaj /(box[j].lb())) < prec_bound))	{
					inf_bound[j]=1;
				}
				if (inf_bound[j]==0 && deltaj < delta) 	{
					nexti=j; infnexti=0;delta=deltaj; found =true;
				}
			}

			if (sup_bound[j]==0) {
				deltaj = fabs (primal_solution[j]- box[j].ub());


				if ((fabs (box[j].ub()) < 1 && deltaj < prec_bound) 	||
					(fabs (box[j].ub()) >= 1 && fabs (deltaj/(box[j].ub())) < prec_bound)) {
					sup_bound[j]=1;
				}
				if (sup_bound[j]==0 && deltaj < delta) {
					nexti=j; infnexti=1;delta=deltaj;  found =true;
				}

			}


		}
	}
	else {
		// Default if the primal solution is not available.
		for (int j=0;j<sys.nb_var;j++) {
			if (inf_bound[j]==0) {
				nexti=j;   infnexti=0; found = true;
				break;
			}

			else if  (sup_bound[j]==0) {
				nexti=j;  infnexti=1; found = true;
				break;
			}
		}
	}
	return found;

}


}// end namespace ibex

