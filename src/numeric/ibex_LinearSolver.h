/*
 * ibex_LinearSolver.h
 *
 *  Created on: 15 mai 2013
 *      Author: nininjo
 */

#ifndef IBEX_LINEARSOLVER_H_
#define IBEX_LINEARSOLVER_H_



#define _IBEX_WITH_SOPLEX_ 1

#ifdef _IBEX_WITH_SOPLEX_
#include "soplex.h"
#endif


namespace ibex {

class LinearSolver {

public:

	static const double default_eps;
	static const double default_max_bound;

	typedef enum  {OPTIMAL, INFEASIBLE, UNKNOWN, TIME_OUT, MAX_ITER } Status_Sol;

	typedef enum  {MINIMIZE, MAXIMIZE} Sense;

	typedef enum {OK, FAIL} Status;

	LinearSolver(int nb_vars, int nb_ctr, int max_iter, int max_time_out);

	~LinearSolver();


	Status_Sol solve();

	Status writeFile(std::string name);


// GET
	int getNbRows() const;

	double getObjValue() const;

	Status getCoefConstraint(Matrix& A);

	Status getB(IntervalVector& B);

	Status getPrimalSol(double * prim);

	Status getDualSol(double * dual);

	Status getInfeasibleDir(double * sol);

// SET

	Status cleanConst();

	Status cleanAll();

	Status setMaxIter(int max);

	Status setMaxTimeOut(int time);

	Status setSense(Sense s);

	Status setVarObj(int var, double coef);

	Status initBoundVar(IntervalVector bounds);

	Status setBoundVar(int var, Interval bound);

	Status setEpsilon(double eps);

	Status addConstraint(Vector & row, CmpOp sign, double rhs );


private:

	int nb_ctrs;

	int nb_vars;
	int nb_rows;

	double obj_value;

	double epsilon;

#ifdef _IBEX_WITH_SOPLEX_
	soplex::SoPlex mysoplex;
#endif

};


} // end namespace ibex

#endif /* IBEX_LINEARSOLVER_H_ */
