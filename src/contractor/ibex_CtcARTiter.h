//============================================================================
//                                  I B E X                                   
// File        : ibex_CtcXNewtonIter.h
// Author      : Ignacio Araya,
//               Bertrand Neveu, Gilles Trombettoni
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jul 20, 2012
// Last Update : March 19, 2013
//============================================================================


#ifndef __IBEX_CTC_ARTITER_H__
#define __IBEX_CTC_ARTITER_H__

#include "ibex_Ctc.h"
#include "ibex_System.h"
#include "ibex_NumConstraint.h"
#include "ibex_CtcLinearRelaxation.h"
#include "ibex_LinearSolver.h"

#include <vector>

namespace ibex {

/** \ingroup ctcgroup
 * \brief ART contractor
 *
 * This class is an implementation of the X-Newton algorithm
 * \author Jordan Ninin
 * \date May 2013
 */

class CtcARTiter : public CtcLinearRelaxation {

public:

	// TODO
	CtcARTiter ();

	~CtcARTiter ();
	// TODO

	/** Basic iteration of the LR-based contractor. Linearize the system and performs calls to Simplex *\
  Apply contraction. It must be implemented in the subclasses **/
	void contract( IntervalVector& box);

	/** ART iteration.
  Linearize the system and performs 2n calls to Simplex in order to reduce
  the 2 bounds of each variable */
	void Linearization( IntervalVector & box);


};

} // end namespace ibex


#endif /* CTC_ARTITER_H_ */

