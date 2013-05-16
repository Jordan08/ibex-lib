//============================================================================
//                                  I B E X                                   
// File        : ART : fixpoint of compo of ARTIter and contractor
// Author      : Jordan Ninin
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 15, 2013
// Last Update : May 15, 2013
//============================================================================

#include "ibex_CtcXNewton.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_CtcCompo.h"
using namespace std;
namespace ibex {

/*! Default fixpoint ratio. */
  const double CtcXNewton::default_xnewton_ratio = 0.2;

  CtcXNewton::CtcART (CtcARTIter & artiter, Ctc & ctc, double ratio) :
    CtcFixPoint(* (new CtcCompo (artiter, ctc)), ratio) {}

} // end namespace ibex

