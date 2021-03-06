/* ============================================================================
 * I B E X - Functions
 * ============================================================================
 * Copyright   : Ecole des Mines de Nantes (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Gilles Chabert
 * Created     : Dec 11, 2013
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_FNC_H__
#define __IBEX_FNC_H__

#include "ibex_IntervalMatrix.h"

namespace ibex {

/**
 * \defgroup function Functions
 */

/**
 * \ingroup function
 * \brief Abstract class of function
 */
class Fnc {
public:

	/**
	 * \brief Build the function from R^n to R^m.
	 */
	Fnc(int n, int m);

	/**
	 * \brief Delete this.
	 */
	virtual ~Fnc();

	/**
	 * \brief Return the total number of variables (n).
	 */
	int nb_var() const;

	/**
	 * \brief Return the number of components of f (m).
	 *
	 */
	int image_dim() const;

	/**
	 * \brief Calculate f(box) using interval arithmetic.
	 *
	 * \pre f must be real-valued
	 */
	virtual Interval eval(const IntervalVector& box) const;

	/**
	 * \brief Calculate f(box) using interval arithmetic.
	 */
	virtual IntervalVector eval_vector(const IntervalVector& box) const;

	/**
	 * \brief Calculate f(x) using interval arithmetic.
	 *
	 * \pre f must be matrix-valued
	 */
	virtual IntervalMatrix eval_matrix(const IntervalVector& x) const;

	/**
	 * \brief Calculate the gradient of f.
	 *
	 * \param x - the input box
	 * \param g - where the gradient has to be stored (output parameter).
	 *
	 * \pre f must be real-valued
	 */
	virtual void gradient(const IntervalVector& x, IntervalVector& g) const;

	/**
	 * \brief Calculate the gradient of f.
	 * \pre f must be real-valued
	 */
	IntervalVector gradient(const IntervalVector& x) const;

	/**
	 * \brief Calculate the Jacobian matrix of f
	 *
	 * \param x - the input box
	 * \param J - where the Jacobian matrix has to be stored (output parameter).
	 *
	 */
	virtual void jacobian(const IntervalVector& x, IntervalMatrix& J) const;

	/**
	 * \brief Calculate the Jacobian matrix of f
	 * \pre f must be vector-valued
	 */
	IntervalMatrix jacobian(const IntervalVector& x) const;

	/**
	 * \brief Calculate the Hansen matrix of f
	 */
	void hansen_matrix(const IntervalVector& x, IntervalMatrix& h) const;

	/**
	 * \brief Return the number of used variables
	 */
	int nb_used_vars() const;

	/**
	 * \brief Return the ith used variable
	 *
	 * \pre 0<=i<nb_used_vars().
	 */
	int used_var(int i) const;

protected:
	friend std::ostream& operator<<(std::ostream& os, const Fnc& f);

	/**
	 * \brief Initialize _nb_used_vars and _used_var
	 */
	virtual void generate_used_vars() const;

	/**
	 * \brief Print the function "x->f(x)" (including arguments)
	 */
	virtual void print(std::ostream& os) const;

	/**
	 * \brief Print the expression "f(x)"
	 */
	virtual void print_expr(std::ostream& os) const;

	Fnc(); // temporary construction
	const int _nb_var;
	const int _image_dim;

	// number of used vars (value "-1" means "not yet generated")
	mutable int _nb_used_vars;

public: // TMP
	// only generated if required
	mutable int* _used_var;

};

/**
 * \brief Streams out a function
 */
std::ostream& operator<<(std::ostream& os, const Fnc& f);

/*================================== inline implementations ========================================*/

inline Fnc::Fnc() : _nb_var(0), _image_dim(0), _nb_used_vars(-1), _used_var(NULL) {

}

inline Fnc::Fnc(int n, int m) : _nb_var(n), _image_dim(m), _nb_used_vars(-1), _used_var(NULL) {

}

inline int Fnc::nb_var() const {
	return _nb_var;
}

inline int Fnc::image_dim() const {
	return _image_dim;
}

inline IntervalVector Fnc::gradient(const IntervalVector& x) const {
	IntervalVector g(x.size());
	gradient(x,g);
	return g;
}

inline IntervalMatrix Fnc::jacobian(const IntervalVector& x) const {
	IntervalMatrix J(image_dim(),x.size());
	jacobian(x,J);
	return J;
}

inline int Fnc::nb_used_vars() const {
	if (_nb_used_vars==-1) generate_used_vars();
	return _nb_used_vars;
}

inline int Fnc::used_var(int i) const {
	if (_nb_used_vars==-1) generate_used_vars();
	return _used_var[i];
}

inline std::ostream& operator<<(std::ostream& os, const Fnc& f) {
	f.print(os);
	return os;
}


} // end namespace ibex

#endif // __IBEX_FNC_H__
