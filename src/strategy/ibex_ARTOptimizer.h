//============================================================================
//                                  I B E X                                   
// File        : ibex_ARTOptimizer.h
// Author      : Jordan Ninin
// License     : See the LICENSE file
// Created     : May 20, 2013
// Last Update : May 20, 2013
//============================================================================

#ifndef __IBEX_ART_OPTIMIZER_H__
#define __IBEX_ART_OPTIMIZER_H__

#include "ibex_Optimizer.h"
#include "ibex_CtcCompo.h"
#include "ibex_CtcART.h"
#include "ibex_CtcARTiter.h"

namespace ibex {

/**
 * \ingroup strategy
 * \brief Default optimizer.
 */
class ARTOptimizer : public Optimizer {
public:
	/**
	 * \brief Create a default optimizer.
	 *
	 * \param sys       - The system to optimize
	 * \param prec      - Stopping criterion for box splitting (absolute precision)
	 * \param goal_prec - Stopping criterion for the objective (relative precision)
	 */
    ARTOptimizer(System& sys, double prec, double goal_prec);

	/**
	 * \brief Delete *this.
	 */
    ~ARTOptimizer();

private:
	Array<Ctc>*  contractor_list (System& sys, System& ext_sys,double prec);


	// -------- information stored for cleanup ----------
	// Extended system
	// (the objective is added in the system as the last variable and the first constraint
    // is used for contraction and bisection)
    System* __ext_sys;

	CtcCompo* __ctc;
	Bsc* __bsc;
};

} // end namespace ibex

#endif // __IBEX_ART_OPTIMIZER_H__
