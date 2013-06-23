//============================================================================
//                                  I B E X                                   
// File        : ibex_Sweep.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jan 22, 2013
// Last Update : Jan 22, 2013
//============================================================================

#include <cassert>
#include "ibex_Sweep.h"

namespace ibex {

Sweep::Sweep(const System& csp, int* order, double jump_ratio) : Ctc(csp.nb_var), order(order), csp(csp),
		jump_ratio(jump_ratio), jump_vector(csp.nb_var), _nb_jumps(0), _nb_trials(0), wl(csp.nb_ctr) {
	assert(csp.nb_var>1);
}


void Sweep::contract(IntervalVector& box) {
	bool res = sweep(box,true);
	res |= sweep(box,false);
	if (res) cout << "pruning successful!\n";
}

//
// We can inactivate the waking list strategy and
// use a simple round-robin strategy instead (this
// strategy is proposed in comments).
//
bool Sweep::try_inflate(const Vector& pt, IntervalVector& X, const Vector& min_width) {

	//int m=csp.nb_ctr();  // for [basic round-robin]
	int c=ctr_num;

	do {
		//cout << "\n\ntry inflation with ctr " << c << ": " << csp.ctrs[c] << "\n";
		//cout << "   pt=" << pt << endl;
		//cout << "   dom=" << X << endl;
		_nb_trials++;

		// calculate the image of the point
		Interval y=csp.f[c].eval(pt);

		// find the interval corresponding to the "negation" of the constraint
		Interval neg;
		switch (csp.ctrs[c].op) {
		case LT:  neg=Interval::POS_REALS; break;
		case LEQ: neg=Interval(next_float(0),POS_INFINITY); break;
		case EQ:
			if (y.lb()<0)
				neg=Interval(next_float(0),POS_INFINITY);
			else
				neg=Interval(NEG_INFINITY,previous_float(0));
			break;
		case GEQ: neg=Interval(NEG_INFINITY,previous_float(0)); break;
		case GT:  neg=Interval::NEG_REALS; break;
		}

		// see if pt is inside the forbidden region of c
		if (y.is_subset(neg)) {

			// try to inflate with c
			csp.f[c].iproj(neg,X,pt);

			bool success = true;
			for (int i=0; i<X.size(); i++)
				success &= (X[i].diam() >= min_width[i]);

			// if the expansion is significant enough, return the box
			if (success) {
				//cout << c << "   ok!" << wl << endl;
				_nb_jumps++;       // only for statistics
				//ctr_num = (c+1)%m; // [basic round-robin] update constraint number for next call
				return true;
			} else {
				//cout << c << "   insignificant " << wl << endl;
			}
		}
		else
			//cout << "   feasible (abort)\n"
			;
		//c=(c+1)%m;               // [[basic round-robin]]
		c = wl.next_candidate(WakingList::NO_JUMP);
		//  } while (c!=ctr_num); // [[basic round-robin]]
	} while (c!=-1);

	// do not replace this "do...while" by a "for" loop:
	// whatever is the initial value of ctr_num, the condition c!=ctr_num would
	// prevent the corresponding constraint to be treated the first time try_inflate is called.

	// at this point, X is pt (possibly slightly enlarged)
	//cout << "\n  no more expansion...\n";
	return false;
}


bool Sweep::sweep(IntervalVector& box, bool min) {

	//  ctr_num = 0;     // [[basic round-robin]]
	ctr_num = wl.first_candidate();

	int k=csp.nb_var;

	bool b=true;         // feasibility?
	bool success=false;  // contraction?

	// The "working" area is the box inside which we try to build forbidden boxes.
	// This concept allows to calculate the actual rate of improvement in try_inflate
	// (since a wide box inside [box] may actually be very thin inside [area]).
	//
	// Note that the upper right corner of the working area has a different meaning
	// than in the discrete version of sweeping: it is not the potentially lowest
	// feasible vector but the potentially highest infeasible one. */

	IntervalVector area(box);

	Vector min_width = jump_ratio * box.diam();

	//cout << "============================ prune_min =============================" << endl;

	while (b && try_inflate(min? area.lb() : area.ub(), area, min_width)) {

		//cout << "  sweep point=" << area.lb() << endl;
		//cout << "  area=" << area << endl;

		b=false;
		bool main_jump=false;
		for (int j=k-1; j>=0; j--) {
			int j2 = order[j]; // simpler alternative: (j+d) % k;

			if ((min  && area[j2].ub() < box[j2].ub()) ||
			    (!min && area[j2].lb() > box[j2].lb())) {
				if (j==0) {
					main_jump = true;
					success = true;
					//ctr_num = rand() % csp.nb_ctr; // [[basic round-robin]] introduce diversity when jumping on main dimension
				}
				area[j2] = min ? Interval(area[j2].ub(), box[j2].ub()) :
							     Interval(box[j2].lb(),  area[j2].lb());
				b = true;
				break;
			} else {
				//cout << "dimension " << j2 << " saturated...\n";
				area[j2] = box[j2];
			}
		}
		ctr_num = wl.next_candidate(main_jump? WakingList::MAIN_JUMP : WakingList::JUMP);
		if (ctr_num==-1) return false;
	}

	if (b) {
		jump_vector = min? area.lb() : area.ub();
		box[order[0]]=Interval(area[order[0]].lb(),box[order[0]].ub());
		return success;
	} else {
		box.set_empty();
		throw EmptyBoxException();
	}
}

} // end namespace
