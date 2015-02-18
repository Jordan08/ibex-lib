//                                  I B E X                                   
// File        : ibex_OptimCtc.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : December 24, 2012
//============================================================================

#include "ibex_OptimCtc.h"
#include "ibex_EmptyBoxException.h"
#include "ibex_Timer.h"

#include "ibex_NoBisectableVariableException.h"

#include <stdlib.h>

using namespace std;

namespace ibex {

const double OptimCtc::default_prec = 1e-07;
const double OptimCtc::default_goal_rel_prec = 1e-07;
const double OptimCtc::default_goal_abs_prec = 1e-07;



OptimCtc::OptimCtc( Ctc& ctc_out, Ctc&  ctc_in, Function& f_cost, Bsc& bsc, double prec,
		double goal_rel_prec, double goal_abs_prec) :
			n(f_cost.nb_var()),
			_ctc_out(ctc_out), _ctc_in(ctc_in), _f_cost(f_cost),
			_localopti(f_cost, IntervalVector(f_cost.nb_var())),
			bsc(bsc), buffer(n),  // first buffer with LB, second buffer with crit (default UB))
			prec(prec), goal_rel_prec(goal_rel_prec), goal_abs_prec(goal_abs_prec),
			trace(false), timeout(1e08), time(0), nb_cells(0),
			loup(POS_INFINITY), uplo(NEG_INFINITY),	loup_point(n),
			loup_changed(false),  uplo_of_epsboxes(POS_INFINITY) {


	// =============================================================

	if (trace) cout.precision(12);
}

OptimCtc::~OptimCtc() {

	buffer.flush();

}

// compute the value ymax (decreasing the loup with the precision)
// the heap and the current box are contracted with y <= ymax
double OptimCtc::compute_ymax() {
	double ymax= loup - goal_rel_prec*fabs(loup);
	if (loup - goal_abs_prec < ymax)
		ymax = loup - goal_abs_prec;
	return ymax;
}


// 2 methods for searching a better feasible point and a better loup

bool OptimCtc::localsearch(const IntervalVector* box, int nb) {
	if (nb<=0) return false;

	bool loup_change=false;
	{
		Vector v_tmp(box[0].size());
		for (int i=0; i<nb;i++) {

			_localopti.set_box(box[i]);
			UnconstrainedLocalSearch::ReturnCode code =
					_localopti.minimize(box[i].random(),v_tmp,goal_rel_prec,100);
			//cout << "[localsearch] "  <<code << " box " << box[i]   << endl;
			if (code == UnconstrainedLocalSearch::SUCCESS) {
				Interval tmp = _f_cost.eval(v_tmp);
				if (tmp.ub()<loup) {
					//update the loup
					loup = tmp.ub();
					loup_point = v_tmp;
					loup_change = true;
					{
						int prec=cout.precision();
						cout.precision(12);
						cout << "[localsearch]"  << " loup update " << loup  << " loup point  " << loup_point << endl;
						cout.precision(prec);
					}
				}
			}
		}
	}
	return loup_change;

}

// a feasible point is a point contracted by ctc_in
bool OptimCtc::direct_try( const Vector point) {
	bool loup_change=false;
	Interval tmp = _f_cost.eval(point);
	if (tmp.ub()<loup) {
		BitSet flags(BitSet::empty(Ctc::NB_OUTPUT_FLAGS));
		try {
			IntervalVector tmpv(point);
			_ctc_out.contract(tmpv, BitSet::all(_ctc_in.nb_var),flags);  //  <---
		} catch (EmptyBoxException &) { 	}  //  <---

		if (flags[Ctc::INACTIVE]) {
			//update the loup
			loup = tmp.ub();
			loup_point = point;
			loup_change = true;

			int prec=cout.precision();
			cout.precision(12);
			cout << "[direct out]"  << " loup update " << loup  << " loup point  " << loup_point << endl;
			cout.precision(prec);
		} else {
			try {
				IntervalVector tmpv(point);
				_ctc_in.contract(tmpv);  //  <---

			} catch (EmptyBoxException &) {
				//update the loup
				loup = tmp.ub();
				loup_point = point;
				loup_change = true;

				int prec=cout.precision();
				cout.precision(12);
				cout << "[direct in]"  << " loup update " << loup  << " loup point  " << loup_point << endl;
				cout.precision(prec);


			}  //  <---
		}
	}
	return loup_change;
}



double minimum (double a, double b) {
	if(a<=b) return a;
	else return b;
}


void OptimCtc::update_uplo() {
	double new_uplo=POS_INFINITY;

	if (! buffer.empty()){
		new_uplo= buffer.minimum();
		if (new_uplo < uplo_of_epsboxes) uplo = new_uplo;
		else uplo= uplo_of_epsboxes;
	}
	else if (buffer.empty() && loup != POS_INFINITY) {
		// empty buffer : new uplo is set to ymax (loup - precision) if a loup has been found
		new_uplo=compute_ymax(); // not new_uplo=loup, because constraint y <= ymax was enforced
		//    cout << " new uplo buffer empty " << new_uplo << " uplo " << uplo << endl;

		double m = minimum(new_uplo, uplo_of_epsboxes);
		if (uplo < m) uplo = m; // warning: hides the field "m" of the class
		// note: we always have uplo <= uplo_of_epsboxes but we may have uplo > new_uplo, because
		// ymax is strictly lower than the loup.
	}

}



/* contract the box of the cell c , try to find a new loup :;
     push the cell  in the 2 heaps or if the contraction makes the box empty, delete the cell. For diversification, rebuild the 2 heaps
 */

void OptimCtc::handle_cell(OptimCell& c, const IntervalVector& init_box ){
	try{
		compute_pf(c);
		contract_and_bound(c, init_box);  // may throw EmptyBoxException

		if (loup < 1.e8)
			c.loup=loup;
		else
			c.loup=1.e8;

		// the cell is put into the 2 heaps
		buffer.push_costpf(&c);

		nb_cells++;

	}
	catch(EmptyBoxException&) {
		draw_vibes(c.box,IntervalVector(1,Interval::EMPTY_SET),"r");
		delete &c;
	}
}

void OptimCtc::compute_pf(OptimCell& c) {
	c.pf &=_f_cost.eval(c.box);
	//c.pf &= _f_cost.eval_affine2(c.box);
}


void OptimCtc::add_buffer_pf(IntervalVector* list, int size) {
	OptimCell* cell;
	double ymax;
	Interval tmp;

	if (loup==POS_INFINITY) ymax=POS_INFINITY;
	else ymax= compute_ymax()+1.e-15;

	for (int i=0; i<size;i++) {

		tmp = _f_cost.eval(list[i]);
		try {
			if (ymax <=tmp.ub()) {
				// Contract with  f_cost(list[i]) <= ymax
				tmp &= Interval(NEG_INFINITY,ymax);
				HC4Revise(AFFINE_MODE).proj(_f_cost,Domain(tmp),list[i]); /// equivalent to :_f_cost.backward(tmp,list[i]);
			}
			cell=new OptimCell(list[i]);
			cell->pu  = 1;
			cell->pf  = tmp;
			cell->loup= loup;
			// add data required by the bisector
			bsc.add_backtrackable(*cell);

			// the cell is put into the 2 heaps with the cost stored in pf
			buffer.push_costpf(cell);

		} catch (EmptyBoxException& ) {/* nothing to do */}
	}

}



void OptimCtc::contract_and_bound(OptimCell& c, const IntervalVector& init_box) {

	/*======================== contract y with y<=loup ========================*/

	double ymax;
	if (loup==POS_INFINITY) ymax=POS_INFINITY;
	else ymax= compute_ymax()+1.e-15;


	// Contract with  f_cost(c.box) <= loup
	if (ymax <=c.pf.ub()) {
		try {
			IntervalVector save_box(c.box);
			c.pf &= Interval(NEG_INFINITY,ymax);
			HC4Revise(AFFINE_MODE).proj(_f_cost,Domain(c.pf),c.box); /// equivalent to :_f_cost.backward(c.pf,c.box);
			draw_vibes(save_box,c.box,"[g]");
		} catch (EmptyBoxException& e) {
			draw_vibes(c.box,IntervalVector(1,Interval::EMPTY_SET),"[g]");
			c.box.set_empty();
			throw e;
		}
	}

	/*================ contract x with Ctc_out ================*/
	//cout << " [contract_out]  x before=" << c.box << endl;
	//cout << " [contract_out]  y before=" << y << endl;

	// Contract only if all the constraint are not satisfied
	if (c.pu ==0)	{

		BitSet flags(BitSet::empty(Ctc::NB_OUTPUT_FLAGS));
		_ctc_out.contract(c.box, BitSet::all(_ctc_out.nb_var),flags); // no need to catch, it is done in handle_cell

		if (flags[Ctc::INACTIVE]) {
			// all the constraint are inactive ( = the entire box is feasible)
			//cout << " [contract_out] box entirely feasible " << endl;
			c.pu =1;
		}

	}
	//cout << " [contract_out]  x after=" << c.box << endl;
	//cout << " [contract_out]  y after=" << y << endl;
	/*====================================================================*/


	/*================ contract x with Ctc_in ================*/
	//cout << " [contract_in]  x before=" << c.box << endl;
	//cout << " [contract_in]  y before=" << y << endl;

	// Contract only if all the constraint are not satisfied

	IntervalVector * box_ok;
	int size_box_ok = 0;
	if (c.pu ==0)	{
		try {
			IntervalVector tmp(c.box);
			BitSet flags(BitSet::empty(Ctc::NB_OUTPUT_FLAGS));
			_ctc_in.contract(tmp, BitSet::all(_ctc_in.nb_var),flags);
			if (tmp==c.box) {
				if (flags[Ctc::INACTIVE]) {
				// all the constraint are active ( = the entire box is unfeasible)
					//cout << " [contract_in] box entirely infeasible " << endl;
					c.box.set_empty();
					throw EmptyBoxException();
				}
			} else {
				size_box_ok = c.box.diff(tmp,box_ok);
				if ((size_box_ok==1)&&(box_ok[0].is_empty())) {
					size_box_ok=0;
					delete[] box_ok;
				}
				else {
					c.box = tmp;
				}

			}


		} catch (EmptyBoxException &) {
			//cout << " [contract_in] box entirely feasible " << endl;
			//draw_vibes(c.box,IntervalVector(1,Interval::EMPTY_SET),"b");
			c.pu=1;
		}
	}

	//cout << " [contract_in]  x after=" << c.box << endl;
	//cout << " [contract_in]  y after=" << y << endl;
	/*====================================================================*/

	/*========================= update loup =============================*/

	bool loup_ch = localsearch(box_ok,size_box_ok);
	if (c.pu==1) {
		//cout<<c.box<<endl;
		loup_ch = (loup_ch || localsearch(&(c.box),1));
	}
	else  {
		// try on the middle point
		loup_ch = (loup_ch || direct_try(c.box.random()));
	}

	c.pf &=_f_cost.eval_affine2(c.box);
	// update of the upper bound of y in case of a new loup found
	if (loup_ch)  	c.pf &= Interval(NEG_INFINITY,compute_ymax());
	loup_changed |= loup_ch;

	if (c.pf.is_empty()) {
		c.box.set_empty();
		throw EmptyBoxException();
	}
	/*====================================================================*/

	/*=====================add box_ok in the buffer=======================*/
	if (size_box_ok>0) {
		add_buffer_pf(box_ok,size_box_ok);
		delete [] box_ok;
		size_box_ok=0;
	}
	/*====================================================================*/


	if ((c.box.max_diam()<=prec) || !c.box.is_bisectable()) {
		// rem1: tmp_box and not c.box because y is handled with goal_rel_prec and goal_abs_prec
		// rem2: do not use a precision contractor here since it would make the box empty (and y==(-inf,-inf)!!)
		// rem 3 : the extended  boxes with no bisectable  domains  should be catched for avoiding infinite bisections
		update_uplo_of_epsboxes(c.pf.lb());

		throw EmptyBoxException();
	}

	// monotonicity_analysis
	if (c.pu ==1) {
		IntervalVector g(n);
		_f_cost.gradient(c.box,g);
		for (int j=0; j<n; j++) {
			if      ((g[j].lb()>=0)&&(c.box[j].lb()!=NEG_INFINITY)) c.box[j]=c.box[j].lb();
			else if ((g[j].ub()<=0)&&(c.box[j].ub()!=POS_INFINITY)) c.box[j]=c.box[j].ub();
		}

	}

}



void OptimCtc::contract ( IntervalVector& box) {
	//IntervalVector save_box(box);
	_ctc_out.contract(box);
	//draw_vibes(save_box,box,"[r]");
}

void OptimCtc::optimize(const IntervalVector& init_box, double obj_init_bound) {
	loup=obj_init_bound;
	if (trace >= 2) cout << "--START--"<< endl;
	uplo=NEG_INFINITY;
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	OptimCell* root=new OptimCell(init_box);
	root->pu = 0;
	root->pf = _f_cost.eval_affine2(init_box);
	root->loup=loup;

	// add data required by the bisector
	bsc.add_backtrackable(*root);

	loup_changed=false;
	loup_point=Vector(init_box.size());
	time=0;
	Timer::start();
	handle_cell(*root,init_box);
	int indbuf=0;

	try {
		while (!buffer.empty()) {
			if (trace >= 2) cout << " buffer " << ((CellBuffer&) buffer) << endl;

			update_uplo();
			if (buffer.empty()) {
				break;
			}

			loup_changed=false;
			OptimCell *c;

			c=buffer.top();


			try {
				pair<IntervalVector,IntervalVector> boxes=bsc.bisect(*c);
				pair<OptimCell*,OptimCell*> new_cells=c->bisect(boxes.first,boxes.second);

				// Transfer the information to the sons.
				new_cells.first->pf = c->pf;
				new_cells.second->pf= c->pf;
				new_cells.first->pu = c->pu;
				new_cells.second->pu= c->pu;
				new_cells.first->loup = c->loup;
				new_cells.second->loup= c->loup;


				buffer.pop();

				handle_cell(*new_cells.first, init_box);
				handle_cell(*new_cells.second, init_box);

				if (uplo_of_epsboxes == NEG_INFINITY) {
					cout << " possible infinite minimum " << endl;
					break;
				}
				if (loup_changed) {
					// In case of a new upper bound (loup_changed == true), all the boxes
					// with a lower bound greater than (loup - goal_prec) are removed and deleted.
					// Note: if contraction was before bisection, we could have the problem
					// that the current cell is removed by contract_heap. See comments in
					// older version of the code (before revision 284).
					double ymax= compute_ymax();

					buffer.contract_heap(ymax);

					if (ymax <=NEG_INFINITY) {
						if (trace) cout << " infinite value for the minimum " << endl;
						break;
					}
					if (trace) { //setprecision(12);
						cout  << "ymax=" << ymax << " uplo= " <<  uplo << endl;
					}
				}
				update_uplo();
				time_limit_check();

			}
			catch (NoBisectableVariableException& ) {
				draw_vibes(c->box,IntervalVector(1,Interval::EMPTY_SET),"[y]");
				update_uplo_of_epsboxes (c->pf.lb());

				buffer.pop();

				if (c->heap_present==0) delete c;


			}
		}
	}
	catch (TimeOutException& ) {
		return;
	}

	Timer::stop();
	time+= Timer::VIRTUAL_TIMELAPSE();
}

void OptimCtc::update_uplo_of_epsboxes(double ymin) {

	// the current box cannot be bisected.  ymin is a lower bound of the objective on this box
	// uplo of epsboxes can only go down, but not under uplo : it is an upperbound for uplo,
	//that indicates a lowerbound for the objective in all the small boxes
	// found by the precision criterion
	assert (uplo_of_epsboxes >= uplo);
	assert(ymin >= uplo);
	if (uplo_of_epsboxes > ymin) 	{
		uplo_of_epsboxes = ymin;
		if (trace) {
			// it is hard to prove the feasability of a point. So there a lot of small boxes.
			cout << "uplo_of_epsboxes: " <<  uplo_of_epsboxes << " | uplo: " << uplo << " | loup: " << loup << " |"<< endl;
		}
	}
}



void OptimCtc::report() {

	if (timeout >0 &&  time >=timeout ) {
		cout << "time limit " << timeout << "s. reached " << endl;
	}
	// No solution found and optimization stopped with empty buffer  before the required precision is reached => means infeasible problem
	if (buffer.empty() && uplo_of_epsboxes == POS_INFINITY && loup==POS_INFINITY) {
		cout << " infeasible problem " << endl;
		cout << " cpu time used " << time << "s." << endl;
		cout << " number of cells " << nb_cells << endl;
	}

	else {
		cout << " best bound in: [" << uplo << "," << loup << "]" << endl;

		double rel_prec;

		if (loup==POS_INFINITY)
			rel_prec= POS_INFINITY;
		else
			rel_prec=(loup-uplo)/(fabs (loup))-1.e-15;

		double abs_prec=loup-uplo-1.e-15;

		cout << " Relative precision obtained on objective function: " << rel_prec << " " <<
				(rel_prec <= goal_rel_prec? " [passed]" : " [failed]") << "  " << goal_rel_prec <<  endl;

		cout << " Absolute precision obtained on objective function: " << abs_prec << " " <<
				(abs_prec <= goal_abs_prec? " [passed]" : " [failed]") << "  " << goal_abs_prec << endl;
		if (uplo_of_epsboxes != NEG_INFINITY && uplo_of_epsboxes != POS_INFINITY)
			cout << " precision on variable domains obtained " << prec << " "   << " uplo_of_epsboxes " << uplo_of_epsboxes << endl;
		else if (uplo_of_epsboxes == NEG_INFINITY)
			cout << " small boxes with negative infinity objective :  objective not bound " << endl;
		if (loup==POS_INFINITY)
			cout << " no feasible point found " << endl;
		else
			cout << " best feasible point " << loup_point << endl;


		cout << " cpu time used " << time << "s." << endl;
		cout << " number of cells " << nb_cells << endl;
	}

}

/* minimal report for benchmarking */
void OptimCtc::time_cells_report() {
	if (timeout >0 &&  time >=timeout ) {
		cout << "timeout" << timeout << "  " << uplo << " " << loup << " ";}
	else
		cout << time << " " ;
	cout << nb_cells << endl;
}


void OptimCtc::report_perf() {

	double rel_prec;
	if (loup==POS_INFINITY)
		rel_prec= POS_INFINITY;
	else
		rel_prec=(loup-uplo)/(fabs(loup))-1.e-15;

	double abs_prec=loup-uplo-1.e-15;

	cout << (	((rel_prec <= goal_rel_prec)||
			(abs_prec <= goal_abs_prec)||
			((buffer.empty() && uplo_of_epsboxes == POS_INFINITY && loup==POS_INFINITY))||
			(uplo<-1.e300)
	)? " T & " : " F & " );

	cout << uplo << " & " << loup << " & ";
	cout <<  time << "  "<< endl ;
}

void OptimCtc::time_limit_check () {
	Timer::stop();
	time += Timer::VIRTUAL_TIMELAPSE();
	if (timeout >0 &&  time >=timeout ) throw TimeOutException();
	Timer::start();
}


void OptimCtc::draw_vibes( const IntervalVector& X0, const IntervalVector& X,const string color) {
	if (X.is_empty()) {
//		vibes::drawBox(X0,color);
		return;
	}
	if (X==X0) return;     // nothing to draw.
	IntervalVector* rest;
	int n=X0.diff(X,rest); // calculate the set difference
	for (int i=0; i<n; i++) {     // display the boxes
//		vibes::drawBox(rest[i],color);
	}
	delete[] rest;
	return;
}


} // end namespace ibex
