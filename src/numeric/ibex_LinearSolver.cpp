//============================================================================
//                                  I B E X
// Interface with the linear solver
// File        : ibex_LinearSolver.cpp
// Author      : Jordan Ninin
// License     : See the LICENSE file
// Created     : May 15, 2013
// Last Update : May 15, 2013
//============================================================================

#include "ibex_LinearSolver.h"


namespace ibex {

const double LinearSolver::default_eps = 1e-10;
const double LinearSolver::default_max_bound = 1e20;
const int LinearSolver::default_max_time_out=100;
const int LinearSolver::default_max_iter=100;
const double LinearSolver::default_max_diam_box = 1.e6;

#ifdef _IBEX_WITH_SOPLEX_


LinearSolver::LinearSolver(int nb_vars, int nb_ctr, int max_iter, int max_time_out, double eps) :
			nb_ctrs(nb_ctr), nb_vars(nb_vars), nb_rows(0), obj_value(0.0), epsilon(eps) {


	mysoplex= new soplex::SoPlex();

	mysoplex->changeSense(soplex::SPxLP::MINIMIZE);
	mysoplex->setTerminationIter(max_iter);
	mysoplex->setTerminationTime(max_time_out);
	mysoplex->setDelta(epsilon);

	// initialize the number of variables of the LP
	soplex::DSVector col1(0);
	for (int j=0; j<nb_vars; j++){
		mysoplex->addCol(soplex::LPCol(0.0, col1, soplex::infinity, -( soplex::infinity )));
	}

	// initialize the constraint of the bound of the variable
	soplex::DSVector row1(nb_vars);
	for (int j=0; j<nb_vars; j++){
		row1.add (j,1.0);
		mysoplex->addRow(soplex::LPRow(-soplex::infinity, row1,soplex::infinity));
		row1.clear();
	}

	nb_rows += nb_vars;

}

LinearSolver::~LinearSolver() {
	delete mysoplex;
}

LinearSolver::Status_Sol LinearSolver::solve() {

	soplex::SPxSolver::Status stat = soplex::SPxSolver::UNKNOWN;
	LinearSolver::Status_Sol res= UNKNOWN;

	try{
		stat = mysoplex->solve();
		if (stat==soplex::SPxSolver::OPTIMAL) {
			obj_value = mysoplex->objValue();
			res= OPTIMAL;
		}
		else if (stat==soplex::SPxSolver::ABORT_TIME)
			res = TIME_OUT;
		else if (stat==soplex::SPxSolver::ABORT_ITER)
			res = MAX_ITER;
		else if (stat==soplex::SPxSolver::INFEASIBLE)
			res = INFEASIBLE;
		else
			res = UNKNOWN;

	}catch(soplex::SPxException&){
		res = UNKNOWN;
	}
//std::cout <<"	stat soplex  "<<stat<<"   "<<res << std::endl;
	return res;

}

LinearSolver::Status LinearSolver::writeFile(std::string name) {
	LinearSolver::Status res= FAIL;
	try {
		mysoplex->writeFile("dump.lp", NULL, NULL, NULL);
		res = OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

int LinearSolver::getNbRows() const {
	return nb_rows;
}

double LinearSolver::getObjValue() const {
	return obj_value;
}


double LinearSolver::getEpsilon() const {
	return epsilon;
}

LinearSolver::Status LinearSolver::getCoefConstraint(Matrix &A) {
	LinearSolver::Status res= FAIL;
	try {
		for (int i=0;i<nb_rows; i++){
			for (int j=0;j<nb_vars; j++){
				A.row(i)[j] =  mysoplex->rowVector(i)[j];
			}
		}
		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::getCoefConstraint_trans(Matrix &A_trans) {
	LinearSolver::Status res= FAIL;
	try {
		for (int i=0;i<nb_rows; i++){
			for (int j=0;j<nb_vars; j++){
				A_trans.row(j)[i] =  mysoplex->rowVector(i)[j];
			}
		}
		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status  LinearSolver::getB(IntervalVector& B) {
	LinearSolver::Status res= FAIL;
	try {
		// Get the bounds of the variables
		for (int i=0;i<nb_vars; i++){
			B[i]=Interval( mysoplex->lhs(i) , mysoplex->rhs(i) );
		}

		// Get the bounds of the constraints
		for (int i=nb_vars;i<nb_rows; i++){
			B[i]=Interval( 	(mysoplex->lhs(i)!=-soplex::infinity)? mysoplex->lhs(i):-default_max_bound,
							(mysoplex->rhs(i)!= soplex::infinity)? mysoplex->rhs(i): default_max_bound   );
			//Idea: replace 1e20 (resp. -1e20) by Sup([g_i]) (resp. Inf([g_i])), where [g_i] is an evaluation of the nonlinear function <-- IA
			//           cout << B(i+1) << endl;
		}
		res = OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}


LinearSolver::Status LinearSolver::getPrimalSol(Vector & primal_solution) {
	LinearSolver::Status res= FAIL;
	try {
		// the primal solution : used by choose_next_variable
		soplex::DVector primal(nb_vars);
		mysoplex->getPrimal(primal);
		for (int i=0; i< nb_vars ; i++) {
			primal_solution[i]=primal[i];
		}
		res = OK;
	}
	catch(soplex::SPxException&) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::getDualSol(Vector & dual_solution) {
	LinearSolver::Status res= FAIL;
	try {
		// the dual solution ; used by Neumaier Shcherbina test
		soplex::DVector dual(nb_rows);
		mysoplex->getDual(dual);

		// TODO  WHY ?? Can you justify ?
		for (int i =0; i< nb_rows ; i++)
		{
			if	(	(mysoplex->rhs(i)>=default_max_bound && (dual_solution[i]<0) ) ||
					(mysoplex->lhs(i)<=-default_max_bound && (dual_solution[i]>0) )) {//Modified by IA
				dual_solution[i]=0;
			}
			else
				dual_solution[i]=dual[i];
		}

		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::getInfeasibleDir(Vector & sol) {
	LinearSolver::Status res= FAIL;
	try {
		soplex::SPxSolver::Status stat1;
		soplex::DVector sol_found(nb_rows);
		stat1 = mysoplex->getDualfarkas(sol_found);

		if (stat1==soplex::SPxSolver::OPTIMAL) // TODO I'm not sure of the value that return getDualfarkas
		{

			// TODO  WHY ?? Can you justify ?
			for (int i =0; i< nb_rows ; i++)
			{
				if	(	(mysoplex->rhs(i)>=default_max_bound && (sol_found[i]<0) ) ||
						(mysoplex->lhs(i)<=-default_max_bound && (sol_found[i]>0) )) {//Modified by IA
					sol[i]=0;
				}
				else
					sol[i]=sol_found[i];
			}
			res =OK;
		}
		else
			res = FAIL;

	}
	catch(soplex::SPxException& ) {
		res =  FAIL;
	}

	return res;
}

LinearSolver::Status LinearSolver::cleanConst() {
	LinearSolver::Status res= FAIL;
	try {
		mysoplex->removeRowRange(nb_vars, nb_rows-1);
		nb_rows = nb_vars;
		obj_value = POS_INFINITY;
		res= OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;

}
LinearSolver::Status LinearSolver::cleanAll() {
	LinearSolver::Status res= FAIL;
	try {
		mysoplex->removeRowRange(0, nb_rows-1);
		nb_rows = 0;
		obj_value = POS_INFINITY;
		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}


LinearSolver::Status LinearSolver::setMaxIter(int max) {
	LinearSolver::Status res= FAIL;
	try {
		mysoplex->setTerminationIter(max);
		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::setMaxTimeOut(int time) {
	LinearSolver::Status res= FAIL;
	try {
		double t = time;
		mysoplex->setTerminationTime(t);
		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::setSense(Sense s) {
	LinearSolver::Status res= FAIL;
	try {
		if (s==LinearSolver::MINIMIZE) {
			mysoplex->changeSense(soplex::SPxLP::MINIMIZE);
			res = OK;
		}
		else if (s==LinearSolver::MAXIMIZE) {
			mysoplex->changeSense(soplex::SPxLP::MAXIMIZE);
			res =OK;
		}
		else
			res = FAIL;

	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}



LinearSolver::Status LinearSolver::setVarObj(int var, double coef) {
	LinearSolver::Status res= FAIL;
	try {
		mysoplex->changeObj(var, coef);
		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::initBoundVar(IntervalVector bounds) {
	LinearSolver::Status res= FAIL;
	try {
		for (int j=0; j<nb_vars; j++){
			// Change the LHS and RHS of each constraint associated to the bounds of the variable
			mysoplex->changeRange(j ,bounds[j].lb(),bounds[j].ub());
		}
		res = OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::setBoundVar(int var, Interval bound) {
	LinearSolver::Status res= FAIL;
	try {
		mysoplex->changeRange(var ,bound.lb(),bound.ub());

		//std::cout<< "improve bound var "<<var<< std::endl;
		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::setEpsilon(double eps) {
	LinearSolver::Status res= FAIL;
	try {
		epsilon = eps;
		mysoplex->setDelta(eps);
		res =OK;
	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}
	return res;
}

LinearSolver::Status LinearSolver::addConstraint(ibex::Vector& row, CmpOp sign, double rhs) {
	LinearSolver::Status res= FAIL;
	try {
		soplex::DSVector row1(nb_vars);
		for (int i=0; i< nb_vars ; i++) {
			row1.add(i, row[i]);
		}

		if (sign==LEQ || sign==LT) {
			mysoplex->addRow(soplex::LPRow(-soplex::infinity, row1, rhs));
			nb_rows++;
			res=  OK;
		}
		else if (sign==GEQ || sign==GT) {
			mysoplex->addRow(soplex::LPRow(rhs, row1, soplex::infinity));
			nb_rows++;
			res = OK;
		}
		else
			res = FAIL;

	}
	catch(soplex::SPxException& ) {
		res = FAIL;
	}

	return res;
}







#endif  // END DEF with SOPLEX


/** \brief Stream out \a x. */
std::ostream& operator<<(std::ostream& os, const LinearSolver::Status_Sol x){

	switch(x) {
	case(LinearSolver::OPTIMAL) :{
			os << "OPTIMAL";
			break;
	}
	case(LinearSolver::INFEASIBLE) :{
			os << "INFEASIBLE";
			break;
	}
	case(LinearSolver::TIME_OUT) :{
			os << "TIME_OUT";
			break;
	}
	case(LinearSolver::MAX_ITER) :{
			os << "MAX_ITER";
			break;
	}
	case(LinearSolver::UNKNOWN) :{
		os << "UNKNOWN";
		break;
	}
	}
	return os;

}


} // end namespace ibex


