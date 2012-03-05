/* ============================================================================
 * I B E X - Interval Vector definition
 * ============================================================================
 * Copyright   : Ecole des Mines de Nantes (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Gilles Chabert
 * Created     : Dec 05, 2011
 * ---------------------------------------------------------------------------- */

#include "ibex_IntervalVector.h"
#include <stdlib.h>

namespace ibex {

namespace { // to create anonymous structure/functions

/** \brief Complementary of an Interval
 *
 * Compute the complementary of x. The result is (c1 union c2)
 */
void complI(const Interval& x, Interval& c1, Interval& c2) {
	if (x.is_empty() || x.is_degenerated()) { // x.is_empty() should not happen if called from compl()
		c1=Interval::ALL_REALS;
		c2=Interval::EMPTY_SET;
		return;
	}
	else {
		if (x.lb()>NEG_INFINITY) {
			c1=Interval(NEG_INFINITY,x.lb());
			if (x.ub()<POS_INFINITY)
				c2=Interval(x.ub(),POS_INFINITY);
			else
				c2=Interval::EMPTY_SET;
		} else if (x.ub()<POS_INFINITY) {
			c1=Interval(x.ub(),POS_INFINITY);
			c2=Interval::EMPTY_SET;
		} else {
			c1=c2=Interval::EMPTY_SET;
		}
	}
}

/** \brief x\y
 *
 */
void diffI(const Interval& x, const Interval& y, Interval& c1, Interval& c2) {
	complI(y,c1,c2);
	c1 &= x;
	if (c1.is_degenerated()) c1=Interval::EMPTY_SET;
	c2 &= x;
	if (c2.is_degenerated()) c2=Interval::EMPTY_SET;

	if (c1.is_empty()) {
		c2=c1;
		c2=Interval::EMPTY_SET;
	}
}

} // end namespace


int IntervalVector::diff(const IntervalVector& y, IntervalVector**& result) const {
	const int n=size();
	const IntervalVector& x=*this;
	IntervalVector **tmp = new IntervalVector*[2*n]; // in the worst case, there is 2n boxes
	Interval c1, c2;
	int b=0;
	if (y.is_empty())
		tmp[b++] = new IntervalVector(x); // copy of this
	else
		for (int var=0; var<n; var++) {

			diffI(x[var],y[var],c1,c2);

			if (!c1.is_empty()) {
				IntervalVector& v=*(tmp[b++]=new IntervalVector(n));
				for (int i=0; i<var; i++) v.set(i, y[i]);
				v.set(var, c1);
				for (int i=var+1; i<n; i++) v.set(i, x[i]);

				if (!c2.is_empty()) {
					IntervalVector& v=*(tmp[b++]=new IntervalVector(n));
					for (int i=0; i<var; i++) v.set(i, y[i]);
					v.set(var, c2);
					for (int i=var+1; i<n; i++) v.set(i, x[i]);
				}
			}
		}

	result=new IntervalVector*[b];
	for (int i=0; i<b; i++) result[i]=tmp[i];
	delete[] tmp;

	return b;
}

int IntervalVector::complementary(IntervalVector**& result) const {
	return IntervalVector(size()).diff(*this,result);
}

IntervalVector IntervalVector::random() const {
	IntervalVector b(size());
	for (int i=0; i<size(); i++) {
		const Interval& xi=(*this)[i];
		// get a random number in [-0.5,0.5]
		double r=rand()/(double) RAND_MAX -0.5;
		// pick a point in the domain of the ith variable
		double p = xi.mid() + r*xi.diam();
		// watch dog
		if (p<xi.lb()) p=xi.lb();
		else if (p>xi.ub()) p=xi.ub();
		b.set(i,p);
	}
	return b;
}

} // end namespace