//============================================================================
//                                  I B E X                                   
// File        : ibex_CtcQuantif.cpp
// Author      : Jordan Ninin, Gilles Chabert
// License     : See the LICENSE file
// Created     : Jan 29, 2014
// Last Update : May 7, 2014
//============================================================================

#include "ibex_CtcExist.h"
#include <cassert>

using namespace std;

namespace ibex {

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const ExprSymbol& y5, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4,y5), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const ExprSymbol& y5, const ExprSymbol& y6, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4,y5,y6), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const ExprSymbol& y5, const ExprSymbol& y6, const ExprSymbol& y7, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4,y5,y6,y7), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const ExprSymbol& y5, const ExprSymbol& y6, const ExprSymbol& y7, const ExprSymbol& y8, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4,y5,y6,y7,y8), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const ExprSymbol& y5, const ExprSymbol& y6, const ExprSymbol& y7, const ExprSymbol& y8, const ExprSymbol& y9, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4,y5,y6,y7,y8,y9), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const ExprSymbol& y5, const ExprSymbol& y6, const ExprSymbol& y7, const ExprSymbol& y8, const ExprSymbol& y9, const ExprSymbol& y10, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const ExprSymbol& y5, const ExprSymbol& y6, const ExprSymbol& y7, const ExprSymbol& y8, const ExprSymbol& y9, const ExprSymbol& y10, const ExprSymbol& y11, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11), init_box, prec) {
}

CtcExist::CtcExist(const NumConstraint& ctr,  const ExprSymbol& y1, const ExprSymbol& y2, const ExprSymbol& y3, const ExprSymbol& y4, const ExprSymbol& y5, const ExprSymbol& y6, const ExprSymbol& y7, const ExprSymbol& y8, const ExprSymbol& y9, const ExprSymbol& y10, const ExprSymbol& y11, const ExprSymbol& y12, const IntervalVector& init_box, double prec)
 : CtcQuantif(ctr, Array<const ExprSymbol>(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12), init_box, prec) {
}

CtcExist::CtcExist(Ctc& ctc, const BitSet& vars, const IntervalVector& init_box, double prec, bool own_ctc) :
	CtcQuantif(ctc, vars, init_box, prec, own_ctc) {
}

CtcExist::CtcExist(const NumConstraint& c, const Array<const ExprSymbol>& y, const IntervalVector& y_init, double prec) :
	CtcQuantif(c, y, y_init, prec) {
}

void CtcExist::proceed(const IntervalVector& x_init, const IntervalVector& x_current, IntervalVector& x_res, IntervalVector& y) {
	IntervalVector x = x_current;

	try {
		CtcQuantif::contract(x, y);
	} catch (EmptyBoxException&) {
		return;
	}

	if (!x.is_subset(x_res)) {

		if (y.max_diam()<=prec) {
			x_res |= x;
			if (x_res==x_init) return;
		}
		else {

			l.push(pair<IntervalVector,IntervalVector>(x,y));

			// ============================== sampling =============================
			// To converge faster to the result, we contract with the mid-vector of y.
			// This allows to get an estimate of "res" without waiting for epsilon-sized
			// parameter boxes (getting quickly some estimate is important for pruning).
			try {
				IntervalVector y_mid = y.mid();
				CtcQuantif::contract(x,y_mid);  // x may be contracted here; that's why we pushed it on the stack *before* sampling.
				x_res |= x;
				if (x_res==x_init) return;
			} catch (EmptyBoxException&) {
				// do nothing
			}
			// =======================================================================
		}
	}

}

void CtcExist::contract(IntervalVector& box) {
	assert(box.size()==Ctc::nb_var);

	// the returned box, initially empty
	IntervalVector res=IntervalVector::empty(Ctc::nb_var);

	assert(l.empty()); // even when an exception is thrown by this function, l is empty.

	l.push(pair<IntervalVector,IntervalVector>(box, y_init));
	
	IntervalVector x_save(Ctc::nb_var);
	IntervalVector x(Ctc::nb_var);

	IntervalVector y(nb_param);
	IntervalVector y_mid(nb_param); // for sampling

	while (!l.empty()) {

		// get the domain of variables
		x_save = l.top().first;
		// get and immediately bisect the domain of parameters (strategy inspired by Optimizer)
		pair<IntervalVector,IntervalVector> cut = bsc->bisect(l.top().second);
		
		l.pop();

		// proceed with the two sub-boxes for y
		proceed(box, x_save, res, cut.first);
		proceed(box, x_save, res, cut.second);
	}

	box &= res;
	if (box.is_empty()) throw EmptyBoxException();

}


} // end namespace ibex
