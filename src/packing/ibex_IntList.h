//============================================================================
//                                  I B E X                                   
// File        : ibex_IntList.h
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jan 22, 2013
// Last Update : Jan 22, 2013
//============================================================================

#ifndef IBEX_INTLIST_H_
#define IBEX_INTLIST_H_

#include <iostream>

using std::string;
using std::ostream;

namespace ibex {

/**
 * \ingroup geometry
 *
 * \brief List of integers
 *
 * This class allows to create an manage a list of integers that belong to a fixed range [0..n-1].
 * Each integer can only appear once in the list which makes the list also a unsorted set.
 * The list can be either circular or not. It is implemented with low-level arrays.
 * This is typically useful when a sublist of variables or constraints is required with fast access.
 *
 * Each operation takes O(1), including insertion/removals.
 *
 * The space complexity O(2*n).
 *
 */
class IntList {

public:

	/** \brief Exception related to this class */
	class Exception {
	public:
		virtual void output(ostream&) const=0;
	};

	class NotAnElement : public Exception {
	public:
		NotAnElement(int x) : x(x) {}
		void output(ostream&) const;
		int x;
	};

	class InvalidValue : public Exception {
	public:
		InvalidValue(int x) : x(x) {}
		void output(ostream&) const;
		int x;
	};

	class EmptyList : public Exception {
	public:
		void output(ostream&) const;
	};

	class OutOfBounds : public Exception {
	public:
		void output(ostream&) const;
	};

	class Repetition : public Exception {
	public:
		Repetition(int x) : x(x) {}
		void output(ostream&) const;
		int x;
	};

	/** Create an (initialy empty) list of elements in [0..n-1], either circular or not */
	IntList(int n, bool circular);

	/** Duplicates the list */
	IntList(const IntList&);

	~IntList();

	/** Change the order of the list so that it starts from [new_first]
	 * and ends with the previous element of [new_first].
	 *
	 * \remark In case of a non-circular list the "old" first becomes the next of the "old" last.
	 *
	 * \throw IntListException::INVALID_VALUE if new_first>=n  */
	void reorder(int new_first);

	/** Return the size of the list. */
	int size() const;

	/** Return true if the list is empty */
	bool empty() const;

	/** Return the first integer of the list
	 *
	 *  \throw IntListException::EMPTY_LIST if the list is empty */
	int first() const;

	/** Return the last integer of the list
	 *
	 *  \throw IntListException::EMPTY_LIST if the list is empty */
	int last() const;

	/** Return the element after x in the list (in O(1)).
	 *  Return [first] if the list is circular and x is the last.
	 *
	 * \throw IntListException::INVALID_VALUE if x>=n
	 * \throw IntListException::NOT_AN_ELEMENT if x does not belong to the list
	 * \thorw IntListException::OUT_OF_BOUNDS if the list is not circular and x is the last element */
	int next(int x) const;

	/** Return the element before x in the list (in O(1))
	 *  Return [last] if the list is circular and x is the first.
	 *
	 * \throw IntListException::INVALID_VALUE if x>=n
	 * \throw IntListException::NOT_AN_ELEMENT if x does not belong to the list
	 * \thorw IntListException::OUT_OF_BOUNDS if the list is not circular and x is the first element */
	int prev(int x) const;

	/** Removes an element and returns the following one (or throw IntListException::EMPTY_LIST is the list has got empty).
	 *
	 * \throw IntListException::INVALID_VALUE if x>=n
	 * \throw IntListException::NOT_AN_ELEMENT if x does not belong to the list
	 * \thorw IntListException::OUT_OF_BOUNDS if the list is not circular and x is the last element.
	 * Note that this exception is not raised if the list has got empty!*/
	int remove(int x);

	/** Add an element at the end of the list
	 * \throw IntListException::INVALID_VALUE if x>=n
	 * \throw IntListException::REPETITION if x is already in the list */
	void add_tail(int x);

	/** Add an element in front of the list
	 * \throw IntListException::INVALID_VALUE if x>=n
	 * \throw IntListException::REPETITION if x is already in the list */
	void add_head(int x);

	/** Insert the element y after x in the list
	 *
	 * \throw IntListException::INVALID_VALUE if either x>=n or y>=n
	 * \throw IntListException::REPETITION if y is already in the list
	 * \throw IntListException::NOT_AN_ELEMENT if x does not belong to the list */
	void insert(int x, int y);

	/** Return the maximal size of the list (n) */
	int limit() const { return n; }

	private:
	int n;                // range is [0..n-1]

	int nothing;          // special value used to notify "null pointer" (if next[x]==nothing ==> x is the last of the list)
	int error;            // special value used to notify "does not belong the the list" (if next[x]==error ==> x is not in the list)

	bool circular;        // true if the list is circular

	int _first;           // first element of the list ([nothing] if the list is empty)
	int _last;            // last element of the list ([nothing] if the list is empty)
	int _size;            // size of the list

	int* _next;           // _next[x] points to the following number in the list ([nothing] if x is the last and the list is not circular)
	int* _prev;           // _prev[x] points to the following number in the list ([nothing] if x is the first and the list is not circular)

};

std::ostream& operator<<(std::ostream&, const IntList&);

/** Stream out an integer list exception with details */
std::ostream& operator<<(std::ostream&, const IntList::Exception&);


} // end namespace

#endif
