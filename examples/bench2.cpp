//============================================================================
//                                  I B E X                                   
// File        : bench2.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jan 22, 2013
// Last Update : Jan 22, 2013
//============================================================================

#include <math.h>

#include "ibex.h"
#include "ibex_Sweep.h"

using namespace std;
using namespace ibex;

#define RADIUS      1.0
#define PREC        1e-07
#define JUMP_RATIO  1e-05
#define FREEDOM     IntervalVector(2, Interval(-10*PREC,10*PREC))
#define CRIT_RATIO  0.1

const double ALL_R[]={1.001,2.001,2.156,2.416,2.702,3.001,3.001,3.306,3.614,3.814,
   	   	              3.925,4.145,4.238,4.391,4.551,4.815,4.865,4.866,4.941,5.164};

const int N=18; // the number of circles

const double R=5.1; //ALL_R[N-1];

/*--------------------------------------------------*/

IntervalMatrix sol(N,2);

/*--------------------------------------------------*/

class Sweep2 : public Ctc {
public:

	/**
	 * "sol" must contain the instantiation for
	 * the coordinates of circles n°0.. i-1
	 */
	Sweep2(const IntervalMatrix& sol, Function& dist, int i, bool greedy) : Ctc(2), dist(dist) {
		Variable x(2);
		SystemFactory fac;
		fac.add_var(x);
		for (int j=0; j<i; j++) {
			fac.add_ctr(dist(x,sol[j])>=2*RADIUS);
		}

		IntervalVector zero(2,Interval::ZERO);
		// =========== For greedy heuristic =================
		if (greedy) {
			// 1- place the circle as far as possible from zero
			fac.add_ctr(dist(x,zero)=R-RADIUS);
			// 2- stick this circle to the previous one
			if (i>0) fac.add_ctr(dist(x,sol[i-1])=2*RADIUS);
		} else {
			fac.add_ctr(dist(x,zero)<=R-RADIUS);
		}

		csp = new System(fac);

//		cout << "--------------------------------------------------------\n";
//		cout << "System " << (i+1) << endl;
//		cout << *csp << endl;
//		cout << "--------------------------------------------------------\n";

		hc4 = new CtcHC4(csp->ctrs);
	}

	~Sweep2() {
		delete hc4;
		delete csp;
	}

	void contract(IntervalVector& box) {
		bool loop=true;

		while (loop) {

			loop=false;
			//cout << "before HC4=" << box << endl;
			hc4->contract(box);
			//cout << "after HC4=" << box << endl;
//			{
//				int order[2] = {0,1};
//				Sweep s(*csp, order, JUMP_RATIO);
//				loop |=s.sweep(box,true);
//				loop |=s.sweep(box,false);
//			}
//			//cout << "before HC4=" << box << endl;
//			hc4->contract(box);
//			//cout << "after HC4=" << box << endl;
//			{
//				int order[2] = {1,0};
//				Sweep s(*csp, order, JUMP_RATIO);
//				loop |=s.sweep(box,true);
//				loop |=s.sweep(box,false);
//			}
		}
	}

	Function& dist;
	System *csp;
	CtcHC4 *hc4;
};


/*--------------------------------------------------*/
int main(int argc, char* argv[]) {

	Variable x(2),y(2);
	Function dist(x,y,sqrt(sqr(x[0]-y[0])+sqr(x[1]-y[1])),"dist");

	sol[0][0]=Interval::ZERO;
	sol[0][1]=Interval(R-1,R-1);

	cout << N+1 << endl;                         // print the total number of circles
	cout << 0 << " " << 0 << " " << R << endl;   // print the enclosing circle
	cout << 0 << " " << R-1 << " " << 1 << endl; // print the first circle

	for (int i=1; i<N; i++) {
		bool sol_found=false;

		RoundRobin rr(PREC);
		CellStack buff;
		{
			//==================================================
			// try direct affectation
			//==================================================
			Sweep2 c(sol, dist, i,true);
			IntervalVector box(2,Interval(-R+1,R-1));
			Solver solver(c, rr, buff, PREC);
			// solver.stop_first_solution();
			vector<IntervalVector> sols;
			solver.start(box);

			if (solver.next(sols)) {
				sol_found=true;
				//cout << "circle " << (i+1) << " directly affected!!!\n";
				sol[i] = sols[0] + FREEDOM;
			}
			buff.flush();
		}

		if (!sol_found) {
			//==================================================
			// try minmax algo
			//==================================================
			Sweep2 c(sol, dist, i,false);
			IntervalVector box(2,Interval(-R+1,R-1));
			Solver solver(c, rr, buff, PREC);
			vector<IntervalVector> sols;
			//cout << "initial box=" << box << endl;
			solver.start(box);

			if (solver.next(sols)) {
				sol_found=true;
				//cout << "circle " << (i+1) << " affected\n";
				sol[i] = sols[0] + FREEDOM;
			}
		}

		if (!sol_found) {
			cout << "\n no solution...\n";
			exit(1);
		}

		// prints the circle we have packed
		cout << sol[i][0].mid() << " " << sol[i][1].mid() << " " << 1 << endl;
	}

	// print bounding box
	cout << (-R) << " " << -R << " " << R << " " << R << endl;
}
