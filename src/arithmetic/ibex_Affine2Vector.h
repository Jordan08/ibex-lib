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

#ifndef __IBEX_AFFINE2_VECTOR_H__
#define __IBEX_AFFINE2_VECTOR_H__

#include <cassert>
#include <iostream>
#include <utility>
#include "ibex_Interval.h"
#include "ibex_Affine2.h"
#include "ibex_IntervalVector.h"
#include "ibex_Vector.h"

namespace ibex {

/**
 * \ingroup arithmetic
 *
 * \brief Vector of Affine2 form
 *
 * By convention an empty vector has a dimension. A vector becomes empty
 * when one of its component becomes empty and all the components
 * are set to the empty Interval.
 */
class Affine2Vector {


private:

	Affine2Vector() : _n(0), _vec(NULL) { }

	int _n;             // dimension (size of vec)
	Affine2 *_vec;	   // vector of elements

public:

	/** \brief  Create \a n Affine2 form . All the components are Affine2([-oo,+oo])
	 * \pre n>0
	 */
	explicit Affine2Vector(int n);

	/**
	 * \brief  Create \a n Affine2Vector of dimension \a n with
	 * all the components initialized to \a x.
	 * \pre n>0
	 */
	Affine2Vector(int n, const Interval& x);

	/**
	 * \brief  Create \a n Affine2Vector of dimension \a n with
	 * all the components initialized to \a x.
	 * \pre n>0
	 */
	Affine2Vector(int n, const Affine2& x);

	/**
	 * \brief Create a copy of \a x.
	 */
	Affine2Vector(const Affine2Vector& x);

	/**
	 * \brief Create \a n Affine2Vector  initialized by [bounds[i][0],bounds[i][1]]
	 *
	 * \param bounds an nx2 array of doubles
	 * \pre n>0
	 */
	Affine2Vector(int n, double  bounds[][2]);

	Affine2Vector(const IntervalVector& x);

	/**
	 * \brief Create the degenerated Affine2Vector x
	 *
	 */
	Affine2Vector(const Vector& x);

	/**
	 * \brief Create [empty; ...; empty]
	 *
	 * Create an empty Affine2Vector of dimension \a n
	 * (all the components being empty Intervals)
	 *
	 * \pre n>0
	 */
	static Affine2Vector empty(int n);

	/**
	 * \brief Delete this vector
	 */
	~Affine2Vector();

	/**
	 * \brief Return the ith Affine2
	 *
	 * A return a const reference to the
	 * i^th component (i starts from 0)
	 */
	const Affine2& operator[](int i) const;

	/**
	 * \brief Return the ith Affine2
	 *
	 * A return a non-const reference to the
	 * i^th component (i starts from 0)
	 */
	Affine2& operator[](int i);

	/**
	 * \brief Set this Affine2Vector to the empty Affine2Vector
	 *
	 * The dimension remains the same.
	 */
	void set_empty();

	/**
	 * \brief Set all the elements to 0 (even if empty).
	 *
	 * \note Emptiness is "overridden".
	 */
	void clear();

	/**
	 * \brief Set all the elements to x (even if empty).
	 *
	 * \note Emptiness is "overridden".
	 */
	void init(const Interval& x);

	/**
	 * \brief Add [-rad,+rad] to all the components of *this.
	 *
	 * \return *this.
	 */
	Affine2Vector& inflate(double rad);

	/**
	 * \brief Resize this Affine2Vector.
	 *
	 * If the size is increased, the existing components are not
	 * modified and the new ones are set to (ZERO), even if
	 * (*this) is the empty Interval (however, in this case, the status of
	 * (*this) remains "empty").
	 */
	void resize(int n2);

	/**
	 * \brief Return a subvector.
	 *
	 * \pre (*this) must not be empty
	 * \return [ (*this)[start_index]; ...; (*this)[end_index] ].
	 */
	Affine2Vector subvector(int start_index, int end_index) const;

	/**
	 * \brief Put a subvector into *this at a given position.
	 *
	 * \param start_index - the position where the subvector
	 * \param subvec - the subvector
	 *
	 * \pre (*this) must not be empty
	 */
	void put(int start_index, const Affine2Vector& subvec);

	/**
	 * \brief Assign this Affine2Vector to x.
	 *
	 * \pre Dimensions of this and x must match.
	 * \note Emptiness is overridden.
	 */
	Affine2Vector& operator=(const Affine2Vector& x);

	Affine2Vector& operator=(const IntervalVector& x);

	/**
	 * \brief Return true if the bounds of this Affine2Vector match that of \a x.
	 */
	bool operator==(const Affine2Vector& x) const;


	/**
	 * \brief The dimension (number of components)
	 */
	int size() const;

	/**
	 * \brief Return the IntervalVector compose by the interval of each Affine2 form
	 * \pre (*this) must be nonempty
	 */
	IntervalVector itv() const;


	/**
	 * \brief Return the lower bound vector
	 * \pre (*this) must be nonempty
	 */
	Vector lb() const;

	/**
	 * \brief Return the upper bound vector
	 * \pre (*this) must be nonempty
	 */
	Vector ub() const;

	/**
	 * \brief Return the midpoint
	 * \pre (*this) must be nonempty
	 */
	Vector mid() const;

	/**
	 * \brief Return the mignitude vector.
	 * \pre (*this) must be nonempty
	 */
	Vector mig() const;

	/**
	 * \brief Return the magnitude vector.
	 * \pre (*this) must be nonempty
	 */
	Vector mag() const;

	/**
	 * \brief Return true iff this Affine2Vector is empty
	 */
	bool is_empty() const;

	/**
	 * \brief Return true iff this Affine2Vector is flat.
	 *
	 * An Affine2Vector is "flat" if the radius is 0 on at least one dimension
	 * An empty interval vector is considered as flat.
	 */
	bool is_flat() const;

	/**
	 * \brief True iff this interval vector contains \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of (*this).
	 * \sa #ibex::Interval::contains(double) const.
	 */
	bool contains(const Vector& x) const;

	/**
	 * \brief true iff this interval vector contains an infinite bound.
	 *
	 * \note An empty interval vector is always bounded.
	 */
	bool is_unbounded() const;

	/**
	 * \brief True iff this interval vector is a subset of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.

	 * \note Always return true if this interval vector is empty.

	 * \sa #ibex::Interval::is_subset(const Interval&) const.
	 */
	bool is_subset(const Affine2Vector& x) const;

	/**
	 * \brief True iff this interval vector is a subset of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.

	 * \note Always return true if this interval vector is empty.

	 * \sa #ibex::Interval::is_subset(const Interval&) const.
	 */
	bool is_subset(const IntervalVector& x) const;

	/**
	 * \brief True iff this interval vector is inside the interior of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.
	 *
	 * \note return true if this interval vector is empty and \a x not.
	 *
	 * \sa #ibex::Interval::is_strict_subset(const Interval&) const.
	 */
	bool is_strict_subset(const Affine2Vector& x) const;

	/**
	 * \brief True iff this interval vector is inside the interior of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.
	 *
	 * \note return true if this interval vector is empty and \a x not.
	 *
	 * \sa #ibex::Interval::is_strict_subset(const Interval&) const.
	 */
	bool is_strict_subset(const IntervalVector& x) const;

	/**
	 * \brief True iff this interval vector is a superset of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.

	 * \note Always return true if \a x is empty.

	 * \sa #ibex::Interval::is_superset(const Interval&) const.
	 */
	bool is_superset(const Affine2Vector& x) const;

	/**
	 * \brief True iff this interval vector is a superset of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.

	 * \note Always return true if \a x is empty.

	 * \sa #ibex::Interval::is_superset(const Interval&) const.
	 */
	bool is_superset(const IntervalVector& x) const;

	/**
	 * \brief True iff \a x is inside the interior of (*this).
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.
	 *
	 * \note return true if x is empty and not (*this).
	 *
	 * \sa #ibex::Interval::is_strict_superset(const Interval&) const.
	 */
	bool is_strict_superset(const Affine2Vector& x) const;


	/**
	 * \brief True iff \a x is inside the interior of (*this).
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.
	 *
	 * \note return true if x is empty and not (*this).
	 *
	 * \sa #ibex::Interval::is_strict_superset(const Interval&) const.
	 */
	bool is_strict_superset(const IntervalVector& x) const;


	/**
	 * \brief True iff *this is a vector of zeros.
	 */
	bool is_zero() const;

    /**
     * \brief True iff *this can be bisected along one dimension.
     *
     * \sa #ibex::Interval::is_bisectable().
     */
    bool is_bisectable() const;

    /**
      * \brief Vector of radii.
      */
    Vector rad() const;

    /**
	 * \brief Return the vector of diameters.
	 */
	Vector diam() const;

	/**
	 * \brief Return the index of a component with minimal/maximal diameter.
	 *
	 *  \param min true => minimal diameter
	 *  \throws InvalidAffine2VectorOp if the Affine2Vector is empty.
	 */
	int extr_diam_index(bool min) const;


	/**
	 * \brief Return the indices of all the components, sorted by increasing/decreasing diameter.
	 */
	void sort_indices(bool min, int tab[]) const;

	/**
	 * \brief Return the maximal diameter among all the components.
	 *
	 *  \throws InvalidAffine2VectorOp if the Affine2Vector is empty.
	 */
	double max_diam() const;

	/**
	 * \brief Return the minimal diameter among all the components.
	 *
	 * \throws InvalidAffine2VectorOp if the Affine2Vector is empty.
	 */
	double min_diam() const;

	/**
	 * \brief Return the volume of this interval vector.
	 *
	 * \note Return \c POS_INFINITY if the vector is unbounded and not flat.
	 * \note Return 0 if the vector is flat and not unbounded.
	 * \warning If the interval vector is both flat and unbounded, the result is undefined.
	 * \sa #flat()
	 * \sa #unbounded()
	 */
	double volume() const;

	/**
	 * \brief Return the perimeter of this interval vector.
	 *
	 * \note Return \c POS_INFINITY if unbounded.
	 */
	double perimeter() const;

	 /**
	  * \brief Return max of the delta, for x\subseteq *this [deprecated]
	  *
	  * Deprecated. Kept for compatibility with ibex 1.xx.
	  */
	double maxdelta(const Affine2Vector&);

	/**
	 * \brief Return the relative distance with x.
	 *
	 * \return \f$\displaystyle \max_{i=1..n} rel\_distance([this]_i, x_i)/diam([this]_i)\f$.
	 *
	 * \sa #ibex::distance(const Affine2Vector& x1, const Affine2Vector& x2).
	 * \sa #ibex::Interval::rel_distance(const Interval& x) const.
	 */
	double rel_distance(const Affine2Vector& x) const;


	/**
	 * \brief Return a random vector inside *this.
	 *
	 * \pre (*this) must be nonempty.
	 */
	Vector random() const;

	/**
	 * \brief Bisect the box
	 *
	 * The box is bisected along the dimension \a i
	 * and with a ratio \a ratio. If (*this)[i] is the interval [a,a+d]:
	 * <ul>
	 * <li> The first box of the result is (*this)[0]x...x(*this)[i-1]x[a+ratio*d]x...
	 * <li> The second box is (*this)[0]x...x(*this)[i-1]x[a+ratio*d,a+d]x...
	 * </ul>
	 * Default value for the ratio is 0.5.
	 * \pre 0<ratio<1
	 */
	std::pair<IntervalVector,IntervalVector> bisect(int i, double ratio) const;



};

/** \ingroup arithmetic */
/*@{*/

/**
 * \brief Display the Affine2Vector \a x
 */
std::ostream& operator<<(std::ostream& os, const Affine2Vector& x);

/**
 * \brief Cartesian product of x and y.
 *
 */
Affine2Vector cart_prod(const Affine2Vector& x, const Affine2Vector& y);

/*@}*/

/*============================================ inline implementation ============================================ */

inline Affine2Vector Affine2Vector::empty(int n) {
	return Affine2Vector(n, Interval::EMPTY_SET);
}

inline Affine2Vector::~Affine2Vector() {
	delete[] _vec;
}

inline void Affine2Vector::set_empty() {
	(*this)[0]=Interval::EMPTY_SET;
}

inline const Affine2& Affine2Vector::operator[](int i) const {
	assert(i>=0 && i<_n);
	return _vec[i];
}

inline Affine2& Affine2Vector::operator[](int i) {
	assert(i>=0 && i<_n);
	return _vec[i];
}

inline void Affine2Vector::clear() {
	init(0);
}

inline int Affine2Vector::size() const {
	return _n;
}

inline bool Affine2Vector::is_empty() const {
	return (*this)[0].is_empty();
}

inline bool Affine2Vector::is_superset(const IntervalVector& x) const {
	return x.is_subset((*this).itv());
}

inline bool Affine2Vector::is_superset(const Affine2Vector& x) const {
	return x.is_subset(*this);
}

inline bool Affine2Vector::is_strict_superset(const IntervalVector& x) const {
	return x.is_strict_subset((*this).itv());
}
inline bool Affine2Vector::is_strict_superset(const Affine2Vector& x) const {
	return x.is_strict_subset(*this);
}

inline double Affine2Vector::max_diam() const {
	return (*this)[extr_diam_index(false)].diam();
}

inline double Affine2Vector::min_diam() const {
	return (*this)[extr_diam_index(true)].diam();
}

inline Affine2Vector cart_prod(const Affine2Vector& x, const Affine2Vector& y) {
	Affine2Vector z(x.size()+y.size());
	z.put(0,x);
	z.put(x.size(),y);
	return z;
}

} // end namespace

#endif /* _IBEX_AFFINE2_VECTOR_H_ */