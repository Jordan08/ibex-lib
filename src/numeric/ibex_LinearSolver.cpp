/*
 * ibex_LinearSolver.cpp
 *
 *  Created on: 16 mai 2013
 *      Author: nininjo
 */

#include "ibex_LinearSolver.h"

/*
#ifdef _IBEX_WITH_SOPLEX_
using namespace soplex;
#endif
*/

namespace ibex {

	const double LinearSolver::default_eps = 1e-10;
	const double LinearSolver::default_max_bound = 1e20;


#ifdef _IBEX_WITH_SOPLEX_



LinearSolver::LinearSolver(int nb_vars, int nb_ctr, int max_iter, int max_time_out) :
	nb_ctrs(nb_ctr), nb_vars(nb_vars), nb_rows(0), obj_value(0.0), epsilon(default_eps)  {


	soplex::SoPlex mysoplex= new soplex::SoPlex();

	// initialize the number of variables of the LP
	soplex::DSVector dummycol(0);
	for (int j=0; j<nb_vars; j++){
		mysoplex.addCol(LPCol(0.0, dummycol, soplex::infinity, - soplex::infinity ));
	}

	// initialize the constraint of the bound of the variable
	soplex::DSVector row1(nb_vars);
	for (int j=0; j<nb_vars; j++){
		row1.add (j,1.0);
		mysoplex.addRow(LPRow(-soplex::infinity, row1,soplex::infinity));
		row1.clear();
	}

	nb_rows += nb_vars;

}

LinearSolver::~LinearSolver() {

	delete mysoplex;
}

LinearSolver::Status_Sol LinearSolver::solve() {

	soplex::SPxSolver::Status stat = soplex::SPxSolver::UNKNOWN;

	try{
		stat = mysoplex.solve();
	}catch(soplex::SPxException&){
		stat = soplex::SPxSolver::UNKNOWN;
	}
	
	if (stat==soplex::SPxSolver::OPTIMAL) {
		obj_value = mysoplex.objValue();
		return OPTIMAL;
	}
	else if (stat==soplex::SPxSolver::ABORT_TIME)
		return TIME_OUT;
	else if (stat==soplex::SPxSolver::ABORT_ITER)
		return MAX_ITER;
	else if (stat==soplex::SPxSolver::INFEASIBLE)
		return INFEASIBLE;
	else
		return UNKNOWN;

}

LinearSolver::Status LinearSolver::writeFile(std::string name) {
	mysoplex.writeFile("dump.lp", NULL, NULL, NULL);
	return OK;
}

int LinearSolver::getNbRows() const {
	return nb_rows;
}

double LinearSolver::getObjValue() const {
	return mysoplex.objValue();
}

LinearSolver::Status LinearSolver::getCoefConstraint(Matrix &As) {

	for (int i=0;i<nb_rows; i++){
		for (int j=0;i<nb_vars; i++){
			As.row(i)[j] =  mysoplex.rowVector(i)[j];
		}
	}
	return OK;

}

LinearSolver::Status  LinearSolver::getB(IntervalVector& B) {

	// Get the bounds of the variables
	for (int i=0;i<nb_vars; i++){
		B[i]=Interval( mysoplex.lhs(i) , mysoplex.rhs(i) );
	}

	// Get the bounds of the constraints
	for (int i=nb_vars;i<nb_rows; i++){
		B[i]=Interval( (mysoplex.lhs(i)!=-soplex::infinity)? mysoplex.lhs(i):-default_max_bound,
				       (mysoplex.rhs(i)!= soplex::infinity)? mysoplex.rhs(i): default_max_bound   );
		//Idea: replace 1e20 (resp. -1e20) by Sup([g_i]) (resp. Inf([g_i])), where [g_i] is an evaluation of the nonlinear function <-- IA
		//           cout << B(i+1) << endl;
	}

	return OK;
}


LinearSolver::Status LinearSolver::getPrimalSol(double* primal_solution) {

	try {
		// the primal solution : used by choose_next_variable
		soplex::DVector primal(nb_vars);
		mysoplex.getPrimal(primal);
		for (int i=0; i< nb_vars ; i++) {
			primal_solution[i]=primal[i];
		}
	}
	catch(Exception&) {
		return FAIL;
	}
	return OK;
}

LinearSolver::Status LinearSolver::getDualSol(double* dual_solution) {

	try {
		// the dual solution ; used by Neumaier Shcherbina test
		soplex::DVector dual(nb_rows);
		mysoplex.getDual(dual);

		for (int i=0; i< nb_rows ; i++)
			dual_solution[i]=dual[i];

		// TODO  WHY ?? Can you justify ?
		for (int i =0; i< nb_rows ; i++)
		{
			if	((mysoplex.rhs(i)>=default_max_bound && (dual_solution[i]<0) ) ||
					(mysoplex.lhs(i)<=-default_max_bound && (dual_solution[i]>0) )) {//Modified by IA
				dual_solution[i]=0;
			}
		}
	}
	catch(Exception&) {
		return FAIL;
	}
	return OK;
}

LinearSolver::Status LinearSolver::getInfeasibleDir(double* sol) {
	soplex::SPxSolver::Status stat1;
	soplex::DVector sol_found(nb_rows);
	stat1 = mysoplex.getDualfarkas(sol_found);

	if (stat1==soplex::SPxSolver::OPTIMAL) {
		for (int i=0; i< nb_rows ; i++)
			sol[i]=sol_found[i];

		return OK;
	}
	else
		return FAIL;


}

LinearSolver::Status LinearSolver::cleanConst() {

	mysoplex.removeRowRange(nb_vars, nb_rows-1);
	nb_rows = nb_vars;

	return OK;

}
LinearSolver::Status LinearSolver::cleanAll() {

	mysoplex.removeRowRange(0, nb_rows-1);
	nb_rows = 0;

	return OK;
}


LinearSolver::Status LinearSolver::setMaxIter(int max) {

	mysoplex.setTerminationIter(max);
	return OK;
}

LinearSolver::Status LinearSolver::setMaxTimeOut(int time) {
	double t = time;
	mysoplex.setTerminationTime(t);

	return OK;
}

LinearSolver::Status LinearSolver::setSense(Sense s) {
	if (s==LinearSolver::MINIMIZE)
		mysoplex.changeSense(soplex::SPxLP::MINIMIZE);
	else if (s==LinearSolver::MAXIMIZE)
		mysoplex.changeSense(soplex::SPxLP::MAXIMIZE);
	else
		return FAIL;

	return OK;
}



LinearSolver::Status LinearSolver::setVarObj(int var, double coef) {
	mysoplex.changeObj(var, 1.0);
	return OK;
}

LinearSolver::Status LinearSolver::initBoundVar(IntervalVector bounds) {

	for (int j=0; j<nb_vars; j++){
		mysoplex.changeRange(j ,bounds[j].lb(),bounds[j].ub());
	}

	return OK;
}

LinearSolver::Status LinearSolver::setBoundVar(int var, Interval bound) {

	mysoplex.changeRange(var ,bound.lb(),bound.ub());
	return OK;
}

LinearSolver::Status LinearSolver::setEpsilon(double eps) {
	epsilon = eps;
	mysoplex.setDelta(eps);
	return OK;
}

LinearSolver::Status LinearSolver::addConstraint(ibex::Vector& row, CmpOp sign, double rhs) {


	soplex::DSVector row1(nb_vars);
	if (sign==LEQ || sign==LT) {
		mysoplex.addRow(LPRow(-soplex::infinity, row1, rhs));
		nb_rows++;
		return OK;
	}
	else if (sign==GEQ || sign==GT) {
		mysoplex.addRow(LPRow(rhs, row1, soplex::infinity));
		nb_rows++;
		return OK;
	}
	else
		return FAIL;


}

#endif

} // end namespace ibex


