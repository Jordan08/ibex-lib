//============================================================================
//                                  I B E X                                   
// File        : ibex_Sweep.h
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jan 22, 2013
// Last Update : Jan 22, 2013
//============================================================================

#ifndef __IBEX_SWEEP_H__
#define __IBEX_SWEEP_H__

#include "ibex_System.h"
#include "ibex_Ctc.h"
#include "ibex_CtcHC4.h"

#include "ibex_WakingList.h"

namespace ibex {

/**
 * \ingroup geometry
 *
 * \brief Sweeping contractor.
 */
class Sweep : public Ctc {

public:

	/**
	 * Creates a sweeping contractor for a set of constraints.
	 *
	 * \param csp        - The system
	 * \param order      - Array containing a permutation of the variables fixing the order
	 *                     of sweeping. The variable order[0] is the one to be pruned and
	 *                     order[1],...,order[k] are the other dimensions we sweep over successively.
	 * \param jump_ratio - Criterion for controlling convergence. A forbidden box is considered
	 *                     to be wide enough if, on each dimension, its size is greater than
	 *                     jump_ratio*(the initial interval).
	 */
	Sweep(const System& csp, int* order, double jump_ratio=0.1);

	/**
	 * Applies contraction
	 *
	 * Note: Contraction is only enforced on the "main" variable (order[0]).
	 */
	virtual void contract(IntervalVector& box);

	/**
	 * Prune the lower or upper bound of the "main" component of box (order[0]).
	 *
	 * \param min - True for pruning the lower bound, false for the upper bound.
	 *
	 * \throw #ibex::EmptyBoxException in case of no solution
	 * Return true in case of success (contraction).
	 *
	 * \precondition If min==true, box.lb() must be infeasible.
	 *               If min==false, box.ub() must be infeasible.
	 */
	bool sweep(IntervalVector& box, bool min);

	/**
	 * Must be called after a call to prune_min or prune_max.
	 *
	 * \return The last jump vector found (i.e., the first point that was
	 * not proven to be infeasible).
	 */
	Vector last_jump_vector() { return jump_vector; }

	/**
	 * (only for statistics)
	 *
	 * \return The number of jumps made since this object exists
	 */
	long nb_jumps() { return _nb_jumps; }

	/**
	 * (only for statistics)
	 *
	 * \return The number of jump trials (i.e., the number of
	 * times we have check infeasibility (and potentially called try_inflate()
	 * on a constraint) made since this object exists.
	 *
	 */
	long nb_trials() { return _nb_trials; }

private:

	/**
	 * Try to contract "X" to an (as large as possible) infeasible sub-box containing "pt",
	 * by considering successively the forbidden regions of each constraint in a strategy
	 * encoded in the WakingList "wl". The first constraint chosen by "wl" should be the one
	 * following the last that provoked a successful contraction (i.e., yielded a wide enough
	 * forbidden box), in the previous calls to try_inflate.
	 *
	 * A contraction is considered successful if the diameter of every dimensions of
	 * the contracted box is at least [min_width]. This criterion is absolute here because
	 * try_inflate does not have knowledge of the size of the initial box (a relative ratio
	 * calculation here, such as:
	 *   X_new.rel_distance(X_old) <= 1-jump_ratio ?
	 * could lead to slow convergence, the box X getting smaller with time).
	 *
	 * Once a constraint has led to a significant contraction, the procedure returns.
	 * Hence, only *one* inflation is actually enforced.
	 *
	 * \return true in case of success, false otherwise.
	 */
	bool try_inflate(const Vector& pt, IntervalVector& X, const Vector& min_width);

	/**
	 * The order in which constraints are handled
	 */
	WakingList wl;

	System csp;         // see the corresponding argument of the constructor
	int* order;         // see the corresponding argument of the constructor
	double jump_ratio;  // see the corresponding argument of the constructor

	Vector jump_vector; // see #last_jump_vector()

	/*
	 * Constraints are selected in a round-robin strategy. ctr_num is this number of the last constraint
	 * that has produced a (wide enough) forbidden box.
	 */
	int ctr_num;

	/** for statistics */
	long _nb_jumps;    // see #nb_jumps()
	long _nb_trials;   // see #nb_trials()
};


} // end namespace ibex

#endif // end __IBEX_SWEEP_H__
