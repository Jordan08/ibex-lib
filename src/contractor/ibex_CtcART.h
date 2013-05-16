/*
 * ibex_CtcART.h
 *
 *  Created on: 16 mai 2013
 *      Author: nininjo
 */

#ifndef _IBEX_CTC_ART_H_
#define _IBEX_CTC_ART_H_

#include "ibex_Ctc.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_CtcARTiter.h"



namespace ibex {

/** \ingroup ctcgroup
 * \brief ART contractor
 *
 * This class is an implementation of the fixpoint part of X-Newton algorithm
 *
 *
 */

class CtcART : public CtcFixPoint {

 public:
    /** build a fix of the composition (CtcCompo) of two contractors : an XNewtonIter contractor and another one (typically Hc4)
     *
     * \param xnewtoniter : an ARTiter contractor
     * \param ctc : a contractor
     * \param ratio : the ratio : stopping criterion for the fixpoint (see ibex_CtcFixPoint)
     *
     */
// TODO
    CtcXNewton(CtcARTiter & artiter, Ctc & ctc, double ratio = default_art_ratio);
    ~CtcXNewton() {delete &ctc;};

/** Default ratio used, set to 0.2. */
    static const double default_art_ratio;
  };

}


#endif /* _IBEX_CTC_ART_H_ */
