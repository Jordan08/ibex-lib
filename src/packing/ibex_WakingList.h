//============================================================================
//                                  I B E X                                   
// File        : ibex_WakingList.h
// Author      : Gilles Chabert, Nicolas Beldiceanu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jan 22, 2013
// Last Update : Jan 22, 2013
//============================================================================

#ifndef __IBEX_WAKING_LIST_H__
#define __IBEX_WAKING_LIST_H__

#include "ibex_IntList.h"
#include <iostream>

using namespace std;

namespace ibex {

/**
 * \ingroup geometry
 *
 * \brief Awakening of constraints in a sweep process
 *
 * This class allows to manage a list of constraints such that each call to the next_candidate
 * function selects a constraint according to their past behaviors.
 * Constraints that provided a significant result (JUMP/MAIN_JUMP) are selected
 * in priority w.r.t. to the constraints that provided no result (NO_JUMP).
 *
 * There is actually three types of results a constraint can entail, in increasing order of significance:
 * <ul>
 * <li> NO_JUMP : no result at all
 * <li> JUMP : a jump, but not in the main (outermost) dimension
 * <li> MAIN_JUMP : a jump in the main (outermost) dimension
 * </ul>
 *
 * The constraints that have entailed a JUMP are in an "active list". The others are in a "reserve" list.
 *
 */
class WakingList {

 public:
	/**
	 * \brief Create a waking list for n constraints.
	 */
	WakingList(int nb_ctr);

	/**
	 * \brief Duplicates this waking list (in its current state)
	 */
	WakingList(const WakingList& wl);

	/**
	 * \brief Delete *this.
	 */
	~WakingList();

	typedef enum { NO_JUMP, JUMP, MAIN_JUMP } JumpResult;

	/**
	 * \brief Reinitialize the list and return the first candidate constraint.
	 *
	 * In case of successive full explorations (until next_candidate is -1),
	 * the relative order of constraints in the list is maintained from
	 * one exploration to the other (the active constraints are placed at the beginning
	 * of the reserve).
	 */
	int first_candidate();

	/**
	 * \brief Return the next candidate constraint
	 *
	 * \param res contains the result obtained with the last candidate
	 */
	int next_candidate(JumpResult res);

	/**
	 * \brief Return a reference to the active list
	 */
	const IntList& active_list() const { return active; }

	/**
	 * \brief Return a reference to the reserve list
	 */
	const IntList& reserve_list() const { return reserve; }

 private:

	friend ostream& operator<<(ostream&,const WakingList&);

	int next_active(JumpResult res);
	int next_reserve(JumpResult res);

	int *tag;                      // contains the number of the main jump until which a constraint will remain active
	int main_jump_num;             // counts the number of jumps in the main (outermost) dimension
	bool active_mode;              // states if we are currently exploring the active list (false= reserve list)
	int cur;                       // current constraint (last candidate)

	int loop_counter;
	int last_triggered;

	/*--------------------------------------------------------------------------------
	 */
	/* fields related to the exploration of the active list only (round-robin) */
	// flags for the whole exploration
	int first_triggered;           // first constraint that triggers a jump in the active list (-1 if none)
	// flags for one pass of the exploration
	int pass_triggered_counter;    // number of constraints that have triggered a jump in the current *pass* of the exploration
	int pass_counter;              // number of constraints to be considered in the active list during one pass of exploration
	// (memorized because the size may changed due to removals)
	// one pass must consider once and only once each active constraint

	/* fields related to the exploration of the reserve list only */
	int last_activated;
	/*--------------------------------------------------------------------------------*/

 public:
	IntList active;
	IntList reserve;

};

/**
 * \brief Display the waking list.
 */
ostream& operator<<(ostream& os, const WakingList& w);

} // end namespace ibex

#endif // __IBEX_WAKINGLIST_H__
