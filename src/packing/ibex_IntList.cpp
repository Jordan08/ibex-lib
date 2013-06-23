//============================================================================
//                                  I B E X                                   
// File        : ibex_IntList.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jan 22, 2013
// Last Update : Jan 22, 2013
//============================================================================

#include "ibex_IntList.h"
#include "ibex_Exception.h"

#include <string.h>

namespace ibex {


ostream& operator<<(ostream& os, const IntList::Exception& e) {
	e.output(os);
	return os;
}

void IntList::EmptyList::output(ostream& os) const {
	os << "The list is empty";
}

void IntList::NotAnElement::output(ostream& os) const {
	os << "Not an element :" << x;
}

void IntList::InvalidValue::output(ostream& os) const {
	os << "The number is invalid (out of range):" << x;
}

void IntList::OutOfBounds::output(ostream& os) const {
	os << "Out of bounds";
}

void IntList::Repetition::output(ostream& os) const {
	os << "Element already exists in the list:" << x;
}

ostream& operator<<(ostream& os, const IntList& l) {
	os << "[";
	if (l.size()>0) {
		int cur=l.first();
		while (cur!=l.last()) {
			os << cur << " ";
			cur = l.next(cur);
		}
		os << l.last();
	}
	os << "]";
	return os;
}


/** Create an (initially empty) list of elements in [0..n-1], either circular or not */
IntList::IntList(int n, bool circular) {
	this->n = n;
	this->circular=circular;

	/* we want [error] to be 0xFFFF if an integer is 16 bits, 0xFFFFFFFF if 32 bits, etc. */
	memset(&error, 0xFF, sizeof(int));
	/* we want [nothing] to be error-1 */
	nothing = error-1;

	if (n<0) {
		ibex_error("IntList: the range bound is negative");
	}
	_next=new int[n];
	_prev=new int[n];

	memset(_next, 0xFF, n*sizeof(int)); // nothing in the list initially (i.e., _next[x]=error for all x)
	//  memset(_prev, 0xFF, sizeof(int)*n);  // useless => we never test membership of x with _prev

	_first=nothing;
	_last=nothing;

	_size=0;
}

IntList::IntList(const IntList& l) {
	n=l.n;
	circular=l.circular;
	error=l.error;
	nothing=l.nothing;
	_next=new int[n];
	_prev=new int[n];
	memcpy(_next, l._next, n*sizeof(int));
	memcpy(_prev, l._prev, n*sizeof(int));
	_first=l._first;
	_last=l._last;
	_size=l._size;
}

IntList::~IntList() {
	delete[] _next;
	delete[] _prev;
}


void IntList::reorder(int new_first) {
	if (new_first<0 || new_first>=n) throw InvalidValue(new_first);
	if (_next[new_first]==error) throw NotAnElement(new_first);

	if (new_first==_first) return;

	_next[_last]=_first; // already done if the list is circular
	_prev[_first]=_last; // idem

	_first=new_first;

	_last=_prev[new_first];

	if (!circular) {
		_next[_last]=nothing;
		_prev[_first]=nothing;
	}
}


int IntList::size() const {
	return _size;
}


bool IntList::empty() const {
	return _size==0;
}

int IntList::first() const {
	if (empty()) throw EmptyList();
	return _first;
}

int IntList::last() const {
	if (empty()) throw EmptyList();
	return _last;
}

int IntList::next(int x) const {
	if (x<0 || x>=n) throw InvalidValue(x);
	if (_next[x]==error) throw NotAnElement(x);
	if (_next[x]==nothing) throw OutOfBounds();
	return _next[x];
}

int IntList::prev(int x) const {
	if (x<0 || x>=n) throw InvalidValue(x);
	if (_next[x]==error) throw NotAnElement(x);
	if (_prev[x]==nothing) throw OutOfBounds();
	return _prev[x];
}


int IntList::remove(int x) {
	if (x<0 || x>=n) throw InvalidValue(x);
	if (_next[x]==error) throw NotAnElement(x);

	_size--;                   // there will be one element less

	if (_first==_last) {       // there is only one element
		_first=nothing;        // the sublist is empty
		_last=nothing;         // the sublist is empty
	} else if (x==_first) {    // x is the first element
		_first = _next[x];
		if (circular) {
			_next[_last]=_first;
			_prev[_first]=_last;
		} else {
			_prev[_first]=nothing;
		}
	} else if (x==_last) {    // this was the last element
		_last = _prev[x];
		if (circular) {
			_next[_last]=_first;
			_prev[_first]=_last;
		} else {
			_next[_last]=nothing;
		}
	} else {
		_next[_prev[x]] = _next[x];  // skip x in the list
		_prev[_next[x]] = _prev[x];
	}

	int res=_next[x];              // ignored if the list is empty
	_next[x]=error;
	if (empty()) throw EmptyList();
	if (res==nothing) throw OutOfBounds();
	return res;
}

void IntList::add_tail(int x) {
	if (x<0 || x>=n) throw InvalidValue(x);
	if (_next[x]!=error) throw Repetition(x);

	if (empty()) {         // the sublist was empty
		_first = x;
		_last = x;
		if (circular) {
			_next[x]=x;
			_prev[x]=x;
		}
		else {
			_next[x]=nothing;
			_prev[x]=nothing;
		}
	} else {
		_next[_last]=x;
		_prev[x]=_last;
		if (circular) {
			_next[x]=_first;
			_prev[_first]=x;
		}
		else {
			_next[x]=nothing;
		}
	}

	_size++;

	_last = x;
}

void IntList::add_head(int x) {
	if (x<0 || x>=n) throw InvalidValue(x);
	if (_next[x]!=error) throw Repetition(x);

	if (empty()) {         // the sublist was empty
		_first = x;
		_last = x;
		if (circular) {
			_next[x]=x;
			_prev[x]=x;
		}
		else {
			_next[x]=nothing;
			_prev[x]=nothing;
		}
	} else {
		_next[x]=_first;
		_prev[_first]=x;
		if (circular) {
			_prev[x]=_last;
			_next[_last]=x;
		}
		else {
			_prev[x]=nothing;
		}
	}

	_size++;

	_first = x;
}

void IntList::insert(int x, int y) {

	if (x<0 || x>=n) throw InvalidValue(x);
	if (y<0 || y>=n) throw InvalidValue(y);
	if (_next[x]==error) throw NotAnElement(x);
	if (_next[y]!=error) throw Repetition(y);

	_size++;

	if (x==_last) {
		_next[x]=y;
		_prev[y]=x;
		_last = y;
		if (circular) {
			_next[y]=_first;
			_prev[_first]=y;
		} else {
			_next[y]=nothing;
		}
	} else {
		_next[y]=_next[x];
		_prev[y]=x;

		_prev[_next[x]]=y;
		_next[x]=y;
	}
}


} // end namespace

