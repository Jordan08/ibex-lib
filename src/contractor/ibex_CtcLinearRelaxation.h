//============================================================================
//                                  I B E X
//
// File        : ibex_CtcLinearRelaxation.h
// Author      : Ignacio Araya Bertrand Neveu, Gilles Trombettoni, Jordan Ninin
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jul 20, 2012
// Last Update : May 15, 2013
//============================================================================

#ifndef IBEX_CTCLINEARRELAXATION_H_
#define IBEX_CTCLINEARRELAXATION_H_

#include "ibex_Ctc.h"
#include "ibex_System.h"
#include "ibex_NumConstraint.h"
#include "ibex_LinearSolver.h"

#include <vector>



namespace ibex {

  /** \ingroup ctcgroup
   * \brief Linear relaxation contractor
   *
   * This in an abstract class for linear relaxation-based contractors (AF-based, X-Newton)
   * \author Ignacio Araya
   * \date April 2012
   */

  class CtcLinearRelaxation : public Ctc {

  public:
    typedef enum  {  ONLY_Y, ALL_BOX } ctc_mode;


    /** Creates the CtcLinearRelaxation
     *
     * \param sys The system
     * \param goal_ctr The id of the goal function, -1 in CSPs
     * \param goal   goal function pointer for optimization, NULL for constraint solving
     * \param cmode ALL_BOX (contracts all variables in the box) | ONLY_Y (in optimization :only improves the left bound of the objective)
     * \param max_iter_soplex : the maximum number of iterations for Soplex (default value 100)
     * \param max_diam_box : the maximum diameter of the box for calling Soplex (default value 1.e6)
     */
    CtcLinearRelaxation(const System& sys, int goal_ctr, Function* fgoal,
		  ctc_mode cmode=ALL_BOX, int max_iter=default_max_iter, int time_out=default_max_time_out, double max_diam_box=default_max_diam_box);

    ~CtcLinearRelaxation ()
      { delete mylinearsolver; delete[] primal_solution; if (dual_solution!=NULL) delete[] dual_solution;}

    /** Basic iteration of the LR-based contractor. Linearize the system and performs calls to Simplex *\
    Apply contraction. It must be implemented in the subclasses **/
    virtual void contract( IntervalVector& box);

    /** The system */
    const System& sys;

    /** Default max_diam_box value, set to 1e6  */
    static const double default_max_diam_box;

    /** Default max_time_out, set to 100s  */
    static const int default_max_time_out;

    /** Default max_iter, set to 100 iterations */
    static const int default_max_iter;


    /** The contraint related to the objective function */
    int goal_ctr;

    /** The goal function pointer for optimization, NULL for constraint solving */
    Function* goal;

  protected:

    /** The maximum number of iterations for the linear solver (default value 100 iterations) */
    int max_iter;

    /** The maximum time for the linear solver (default value 100s) */
    int max_time_out;

    /** The maximum diameter of the box for calling the linear solver (default value 1.e6)  */
    double max_diam_box;

    /** Indicates if in optimization  only the objective is contracted (cmode=ONLY_Y) or all the box (ALL_BOX) */
    ctc_mode cmode;

    /** The linear solver that will be use */
    LinearSolver mylinearsolver;

    /* the primal solution found by the LP solver */
    double* primal_solution;

    /* the dual solution found by the LP solver */
    double* dual_solution;

    /** the linearization technique. It must be implemented in the subclasses */
    virtual int Linearization(IntervalVector& box) = 0;

    /*Neumaier Shcherbina postprocessing in case of optimal solution found : the result obj is made reliable */
    void NeumaierShcherbina_postprocessing (int n, int nr, int var, Interval & obj, IntervalVector& box, Matrix & As, IntervalVector& B, bool minimization);

    /* Neumaier Shcherbina postprocessing in case of infeasibilty found by LP  returns true if the infeasibility is proved */
    bool  NeumaierShcherbina_infeasibilitytest (int n, int nr, IntervalVector& box, Matrix & As, IntervalVector& B);

    /* Achterberg heuristic for choosing the next variable  and which bound to optimize */
    void choose_next_variable ( IntervalVector &box,  int & nexti, int & infnexti, int* inf_bound, int* sup_bound);

    /* call to LinearSolver */
    LinearSolver::Status_Sol run_simplex(IntervalVector &box, LinearSolver::Sense sense, int var, int n, Interval & obj, double bound);

    void optimizer(IntervalVector &box, int nb_var, int nb_ctr);

    bool isInner(IntervalVector & box, const System& sys, int j); /* redoundant method? */

  };
}

#endif



#endif /* IBEX_CTCLINEARRELAXATION_H_ */
