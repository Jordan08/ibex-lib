//============================================================================
//                                  I B E X                                   
// File        : ibex_WakingList.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jan 22, 2013
// Last Update : Jan 22, 2013
//============================================================================

#include "ibex_WakingList.h"
#include <assert.h>
#include <string.h>

#define MAX_LOOP 2

namespace ibex {

WakingList::WakingList(int nb_ctr) : tag(new int[nb_ctr]), active(nb_ctr,true), reserve(nb_ctr,false) {
	// all the constraints are initialy in the reserve
	for (int i=0; i<nb_ctr; i++) {
		reserve.add_tail(i);
	}
}

WakingList::WakingList(const WakingList& wl) : active(wl.active), reserve(wl.reserve) {

	int n=active.limit();
	tag = new int[n];
	memcpy(tag, wl.tag, n*sizeof(int));
	main_jump_num = wl.main_jump_num;
	active_mode = wl.active_mode;
	cur = wl.cur;
	first_triggered = wl.first_triggered;
	pass_triggered_counter = wl.pass_triggered_counter;
	pass_counter = wl.pass_counter;
	last_activated = wl.last_activated;
	loop_counter = wl.loop_counter;
	last_triggered = wl.last_triggered;
}

WakingList::~WakingList() {
	delete[] tag;
}

int WakingList::first_candidate() {
	main_jump_num=0;

	// all the constraints of the active list move at
	// the beginning of the reserve
	while (!active.empty()) {
		reserve.add_head(active.last());
		try {
			active.remove(active.last());
		} catch(const IntList::EmptyList&) { }
	}

	active_mode=false;
	last_activated=-1;                        // none yet in the active list
	first_triggered=-1;                       // no main jump yet
	cur = reserve.first();
	pass_triggered_counter=0;                 // in case no constraint is activated in the reserve (avoids loop)
	loop_counter = 0;
	last_triggered = -1;
	return cur;                               // which is 0
}


int WakingList::next_active(JumpResult res) {
	int next;
	pass_counter--;

	if (res==MAIN_JUMP || res==JUMP) {
		tag[cur] = main_jump_num+1;             // will remain active until next sweeping loop

		pass_triggered_counter++;               // another jump occurs in this pass

		if (last_triggered==cur) {
			loop_counter++;
		} else {
			loop_counter=0;
			last_triggered=cur;
		}

		if (first_triggered==-1) {              // a jump occurs in this exploration
			first_triggered = cur;
		}

		if (res==MAIN_JUMP) {
			main_jump_num++;
			pass_counter=0;                       // force restart round-robin
			//pass_triggered_counter=2;           // to avoid going to the reserve
			if (loop_counter>=MAX_LOOP)
				loop_counter=MAX_LOOP-1;          // idem
			next = active.next(first_triggered);  // start from the next of first_triggered (for diversity)
		} else {
			next = active.next(cur);
		}
	}
	else { // NO_JUMP
		if (tag[cur]<main_jump_num) {
			try {
				next=active.remove(cur);
			} catch (const IntList::EmptyList&) {
				// here, we must have pass_counter=0 and first_triggered=-1
				assert(pass_counter==0 && !first_triggered==-1);
			}
			reserve.add_tail(cur);
		} else {
			next=active.next(cur);            // can loop on the same constraint (if there is 1 active constraint)
		}
	}

	if (pass_counter==0) {                    // end of the pass
		if (pass_triggered_counter>0 &&       // the fixpoint is not reached yet
				loop_counter<MAX_LOOP) {      // and we can still accept the same constraint to loop
			pass_triggered_counter = 0;       // reinit for the next pass
			pass_counter = active.size();     // start a new pass
			if (res==MAIN_JUMP) {
				first_triggered = -1;         // reinit the whole exploration
				last_triggered = -1;
				loop_counter= 0;
			}
		} else {
			active_mode = false;
			last_activated = -1;              // init for future insertions from the reserve
			try {
				next=reserve.first();
			} catch (const IntList::EmptyList&) { // the reserve is empty
				return -1;                    // the whole exploration is over
			}
		}
	}

	return next;
}

/** Return the next candidate
 *
 * Return -1 if there is no mode candidate */
int WakingList::next_candidate(JumpResult res) {
	if (active_mode)
		return cur=next_active(res);
	else
		return cur=next_reserve(res);
}

int WakingList::next_reserve(JumpResult res) {

	int next;
	bool end_of_reserve=false;

	if (res==MAIN_JUMP || res==JUMP) {        // if the current constraint has made a jump
		try {
			next = reserve.remove(cur);         // we remove it from the reserve and get the next one in the reserve
		} catch(const IntList::Exception&) {    // the exception is either empty or we have reached the limit
			end_of_reserve=true;
		}
		if (last_activated==-1)                 // and add it in the active list
			active.add_head(cur);               // either at the beginning
		else
			active.insert(last_activated,cur);  // or just after the one previously inserted

		tag[cur]=main_jump_num+1;               // will remain active until next sweeping loop
		if (first_triggered==-1) {              // a jump occurs in this exploration
			first_triggered = cur;
		}

		last_activated=cur;                     // becomes the last activated constraint

		/* the following lines seems to have a bad effect, especially the first time the reserve is explored.
		 * Indeed, in this case, the first constraint in the reserve that performs a jump is
		 * activated, and this constraint is called alone in a loop until it gets inactive. */
		/*    if (res==MAIN_JUMP) {
      active_mode = true;                   // if this was a main jump, we exit the reserve
      return next_active(MAIN_JUMP);        // and restart the exploration of the active list
                                            // everything happens here just as if the constraint
                                            // was already on the main list
					    }
		 */

	} else {
		try {
			next = reserve.next(cur);             // we get the next one in the reserve
		} catch(const IntList::OutOfBounds&) {
			end_of_reserve=true;
		}
	}

	if (end_of_reserve)
		if (last_activated!=-1 ||               // a constraint in the reserve has made a jump
				pass_triggered_counter==1) {    // the active list is not saturated yet
			active_mode = true;                 // we go back to active mode
			pass_triggered_counter = 0;         // reinit for the next pass
			pass_counter = active.size();       // reinit for the next pass
			return active.first();              // start a new pass from the beginning of the active list
		} else {
			return -1;                          // NO JUMP and end of reserve
		}
	else {
		return next;
	}
}

ostream& operator<<(ostream& os, const WakingList& w) {
	os << "<" << w.main_jump_num <<">  active=[";
	const IntList& active=w.active_list();
	if (active.size()>0) {
		int cur=active.first();
		while (cur!=active.last()) {
			os << cur << "(" << w.tag[cur] <<") ";
			cur = active.next(cur);
		}
		cur=active.last();
		os << cur << "(" << w.tag[cur] << ")";
	}
	os << "] reserve=" << w.reserve;
	return os;
}


} // end namespace



