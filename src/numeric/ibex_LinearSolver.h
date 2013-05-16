/*
 * ibex_LinearSolver.h
 *
 *  Created on: 15 mai 2013
 *      Author: nininjo
 */

#ifndef IBEX_LINEARSOLVER_H_
#define IBEX_LINEARSOLVER_H_




namespace ibex {

class LinearSolver {

public:

	typedef enum  {OPTIMAL, INFEASIBLE, UNKNOWN, TIME_OUT, MAX_ITER } Status_Sol;

	typedef enum  {MINIMIZE, MAXIMIZE} Sense;

	typedef enum {OK, FAIL} Status;



	LinearSolver(int nb_vars, int nb_ctr, int max_iter, int max_time_out);

	~LinearSolver();


	Status_Sol solve();

	Status writeFile(std::string name);


// GET
	int getNbRows();

	double getObjValue();

	Matrix& getCoefConstraint();

	IntervalVector& getB();

	double * getPrimalSol();

	double * getDualSol();

	Status getInfeasibleDir(double * sol);

// SET

	Status clean();

	Status cleanObj();

	Status setMaxIter(int max);

	Status setMaxTimeOut(int time);

	Status changeSense(Sense s);

	Status changeVarObj(int var, double coef);

	Status initBoundVar(IntervalVector bounds);

	Status changeBoundVar(int var, Interval bound);

	Status setEpsilon(double eps);

	Status addConstraint(Vector & row, CmpOp sign, double rhs );


private:

	int nb_rows;
	int nb_vars;

	double obj_value;

	double epsilon;


};


} // end namespace ibex

#endif /* IBEX_LINEARSOLVER_H_ */
