//============================================================================
//                                  I B E X                                   
// File        : ibex_ARTOptimizer.cpp
// Author      : Jordan Ninin
// License     : See the LICENSE file
// Created     : May 20, 2013
// Last Update : May 20, 2013
//============================================================================

#include "ibex_ARTOptimizer.h"
#include "ibex_SmearFunction.h"
#include "ibex_CtcHC4.h"
#include "ibex_CtcAcid.h"
#include "ibex_CtcCompo.h"
#include "ibex_BoolMask.h"
#include "ibex_CellStack.h"
#include "ibex_Array.h"
#include <vector>

namespace ibex {

using namespace std;


namespace {

// These variable and function are necessary because we need
// to pass the extended system "ext_sys" to the base class
// constructor before it is built as a member of the class

System* tmp_ext_sys=NULL;

System& get_ext_sys(System& sys) {
	if (tmp_ext_sys==NULL)
		tmp_ext_sys=new System(sys,System::EXTEND);
	return *tmp_ext_sys;
}

}

// the ARTOptimizer constructor  1 point for sample_size
ARTOptimizer::ARTOptimizer(System& _sys, double prec, double goal_prec) :
		Optimizer(_sys,
				  *new SmearSumRelative(get_ext_sys(_sys),prec),
				  *new CtcCompo (* (contractor_list(_sys,get_ext_sys(_sys),prec))), // warning: we don't know which argument is evaluated first (tmp_ext_sys may be NULL)
				  prec,
				  goal_prec,
				  goal_prec,
				  1),
		__ext_sys(tmp_ext_sys),
		__ctc(dynamic_cast<CtcCompo*>(&ctc)), __bsc(&bsc)  {
  
	tmp_ext_sys=NULL; // important! re-init for next call to constructor.

	srand(1);}



// the contractor list  hc4, acid(hc4), xnewton
Array<Ctc>*  ARTOptimizer::contractor_list (System& sys, System& ext_sys,double prec) {
	Array<Ctc>* ctc_list;
	ctc_list= new Array<Ctc>(3);
	// first contractor on ext_sys : incremental hc4  ratio propag 0.01
	ctc_list->set_ref(0, *new CtcHC4 (ext_sys.ctrs,0.01,true));
	// second contractor on ext_sys : acid (hc4)   with incremental hc4  ratio propag 0.1
	ctc_list->set_ref(1, *new CtcAcid (ext_sys,*new CtcHC4 (ext_sys.ctrs,0.1,true),true));
	// the last contractor is CtcXNewtonIter  with rfp=0.2 and rfp2=0.2
	// the limits for calling soplex are the default values 1e6 for the derivatives and 1e6 for the domains : no error found with these bounds
	ctc_list->set_ref(2,*new CtcART
			  (*new CtcARTiter(ext_sys),
			   *new CtcHC4 (ext_sys.ctrs,0.01)));
	return ctc_list;
}


// deletion of all dynamically created objects
ARTOptimizer::~ARTOptimizer() {
	delete &((dynamic_cast<CtcAcid*> (&__ctc->list[1]))->ctc);
	CtcCompo* ctccompo= dynamic_cast<CtcCompo*>(&(dynamic_cast<CtcART*>( &__ctc->list[2])->ctc));

	delete &(ctccompo->list[0]);
	delete &(ctccompo->list[1]);
	for (int i=0 ; i<__ctc->list.size(); i++)
		delete &__ctc->list[i];
	delete __ctc;
	delete __bsc;
	delete __ext_sys;

}

} // end namespace ibex