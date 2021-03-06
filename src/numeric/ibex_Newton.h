//============================================================================
//                                  I B E X                                   
// File        : ibex_Newton.h
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 24, 2012
// Last Update : May 24, 2012
//============================================================================

#ifndef __IBEX_NEWTON_H__
#define __IBEX_NEWTON_H__

#include "ibex_Fnc.h"

namespace ibex {

/**
 * \brief Default Newton precision
 */
extern double default_newton_prec;

/**
 * \brief Default Gauss Seidel ratio
 */
extern double default_gauss_seidel_ratio;

/** \ingroup numeric
 *
 * \brief Multivariate Newton operator (contracting).
 *
 * This function implements the Hansen-Sengupta variant of the contracting interval Newton iteration.
 * It can be used either as a contractor or as an existence test.
 * The underlying linear routine is \link ibex::gauss_seidel(const IntervalMatrix&, const IntervalVector&, IntervalVector&, double) Gauss-Seidel \endlink.
 *
 * \param f - The function
 * \param box - The box
 * \param prec (optional) - Criterion for stopping the main loop of the iteration. If a step of interval Newton does not reduce the
 * variable domain diameter by more than \a prec (this is not ratio but an absolute value), then the procedure stops.
 * The default value is #default_prec (1e-07).
 * \param gauss_seidel_ratio (optional) - Criterion for stopping the inner Gauss-Seidel loop. If a step of Gauss Seidel does not
 * reduce the variable domain diameter by more than \a ratio_gauss_seidel times, then the linear iteration stops.
 * The default value is #default_gauss_seidel_ratio (1e-04).
 * \return True if one variable has been reduced by more than \a prec.
 */
bool newton(const Fnc& f, IntervalVector& box, double prec=default_newton_prec, double gauss_seidel_ratio=default_gauss_seidel_ratio);

/** \ingroup numeric
 *
 * \brief Multivariate Newton operator (inflating).
 *
 * This function implements the inflating Newton iteration as described in page 9 of
 * "Certified Parallelotope Continuation for One-Manifolds" by Martin et al., SINUM 2013.
 *
 * The underlying linear routine is \link ibex::gauss_seidel(const IntervalMatrix&, const IntervalVector&, IntervalVector&, double) Gauss-Seidel \endlink.
 *
 * The box is inflated at each step as follows:
 *    [x] <- mid[x] + delta*(rad[x]) + chi*[-1,+1]
 *
 * \param f                 - The function
 * \param box               - Input/output argument. Input: starting box. Output: In case of success, this box is proven to contain a solution.
 *                            Otherwise, nothing can be said.
 * \param k_max_iteration (optional)
 *                          - maximal number of iterations
 * \param mu_max_divergence (optional)
 *                          - stop criterion. When the Hausdorff distance between two iterates increases by
 *                            a ratio greater than mu_max_divergence, the procedure halts. Value 1.0 by default
 *                            is for detecting divergence.
 * \param delta_relative_inflat (optional)
 *                          - The box is inflated at each step as follows:
 *                            [x] <- mid[x] + delta*(rad[x]) + chi*[-1,+1]
 * \param chi_absolute_inflat (optional)
 *                          - The box is inflated at each step as follows:
 *                            [x] <- mid[x] + delta*(rad[x]) + chi*[-1,+1]
 *
 * \return True if it is proven that the output box contains a solution.
 */
bool inflating_newton(const Fnc& f, IntervalVector& box,
		int k_max_iteration=15, double mu_max_divergence=1.0,
		double delta_relative_inflat=1.1, double chi_absolute_inflat=1e-12);

} // end namespace ibex
#endif // __IBEX_NEWTON_H__
