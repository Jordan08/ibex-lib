//============================================================================
//                                  I B E X
// File        : exo04.cpp
// Author      : Jordan Ninin
// License     : See the LICENSE file
// Created     : Jan 29, 2014
// Last Update : Jan 29, 2014
//============================================================================


#include "ibex.h"
#include <stdlib.h>

using namespace std;
using namespace ibex;

static int ind = 4;

void contract_and_draw( Ctc& c, IntervalVector& X, const string color, const string colorother) {
	IntervalVector X0(X);       // get a copy
	try {
		BitSet flags(BitSet::empty(Ctc::NB_OUTPUT_FLAGS));
		c.contract(X, BitSet::all(c.nb_var),flags);

		if (flags[Ctc::INACTIVE]) {
			vibes::drawBox(X0[ind].lb(),X0[ind].ub(),X0[ind+1].lb(),X0[ind+1].ub(),colorother);
			X.set_empty();
			return;

		}
	} catch(EmptyBoxException&) {
		vibes::drawBox(X0[ind].lb(),X0[ind].ub(),X0[ind+1].lb(),X0[ind+1].ub(),color);
		X.set_empty();
		return;
	}
	if (X==X0) return;     // nothing contracted.
	IntervalVector* rest;
	int n=X0.diff(X,rest); // calculate the set difference
	for (int i=0; i<n; i++) {     // display the boxes
		vibes::drawBox(rest[i][ind].lb(),rest[i][ind].ub(), rest[i][ind+1].lb(),rest[i][ind+1].ub(),color);
	}
	delete[] rest;
	return;
}



int main() {
	{
		std::string figureName = "Avion_0-3";
		vibes::beginDrawing();           // <== Initializes the VIBES "connection"
		vibes::newFigure(figureName);       // <== Create a new VIBes figure

		double epsilon = 5.e-3;
		double gaol_prec= 1.e-2;
		double epsilon_time = 1.e-4;
		double T_final= 1;
		double time_out= 3600;

		Interval proj_time(0,T_final);

		Interval init_time(0.1);
		Interval init_delta(0.4);
		Interval init_q(-0.1,0.1);

		double secu =0.01;
		// number of plane
		int n=2;
		
		// initial position and speed
		Matrix P(n,2); // initial position of each plane
		Matrix V(n,2); // initial speed vector of each plane

		P[0][0] = 0;
		P[0][1] = 0;
		V[0][0] = 1;
		V[0][1] = 1;
/*
		P[1][0] = 1;
		P[1][1] = 0;
		V[1][0] = -1;
		V[1][1] = 1;

		P[2][0] = 0;
		P[2][1] = 1;
		V[2][0] = 1.3;
		V[2][1] = -1;
*/
		P[1][0] = 1;
		P[1][1] = 1;
		V[1][0] = -1;
		V[1][1] = -1.3;


		
		// Symbolic part
		Variable t,tm, delta,p(2),v(2),q, Z(2);
		Variable Tm(n),Delta(n),Q(n);

		Function traj_avion(t,tm,delta,q, p,v,
		  Return(
			chi(t-tm,
				p[0]+v[0]*t,
				chi(t-(tm+delta),
						p[0]-q*v[0]*tm+q*v[0]*t,
						p[0]+q*v[0]*delta + v[0]*t
				  )
			),
			chi(t-tm,
				p[1]+v[1]*t,
				chi(t-(tm+delta),
						p[1]-q*v[1]*tm+q*v[1]*t,
						p[1]+q*v[1]*delta + v[1]*t
				  )
			)
		  )
		);
				                       
		Function norm_2 (Z, sqrt(sqr(Z[0]) +sqr(Z[1])));


		// we want the distance between the robots stay larger than 2.

		vector<Ctc*> array_out;
		for (int i=0; i<n; i++) {
			for(int j=i+1; j<n; j++) {
				NumConstraint *c=new NumConstraint(Tm,Delta,Q,t,
						norm_2(
								traj_avion(t,Tm[i],Delta[i],Q[i],IntervalVector(P[i]),IntervalVector(V[i]))
							   -traj_avion(t,Tm[j],Delta[j],Q[j],IntervalVector(P[j]),IntervalVector(V[j]))
								)
								>=secu);
				array_out.push_back(new CtcFwdBwd(*c));
			}
		}

		// CtcHC4 = (CtcCompo + CtcFixPoint)
		CtcCompo outside1(array_out);
		
		// Create a mask to indicate which variable are projected
		BitSet mask(BitSet::all(3*n+1));
		mask.remove(3*n);  // the last variable "t" will be projected
		// The initial domain of "t" is [0,1]
		CtcForAll outside(outside1,mask,IntervalVector(1,proj_time),epsilon_time);


		// Now we create the contractor for the complementary domaine
		vector<Ctc*> array_in;
		for (int i=0; i<n; i++) {
			for(int j=i+1; j<n; j++) {
				NumConstraint *c= new NumConstraint(Tm,Delta,Q,t,
						norm_2(
								traj_avion(t,Tm[i],Delta[i],Q[i],IntervalVector(P[i]),IntervalVector(V[i]))
							   -traj_avion(t,Tm[j],Delta[j],Q[j],IntervalVector(P[j]),IntervalVector(V[j]))
								)
								<=secu);
				array_in.push_back(new CtcFwdBwd(*c));
			}
		}
		CtcUnion inside1(array_in);
		CtcExist inside(inside1,mask,IntervalVector(1,proj_time),epsilon_time);

		// Build the initial box of V
		IntervalVector init_box(3*n,init_time);
		init_box.put(n,IntervalVector(n,init_delta));
		init_box.put(2*n,IntervalVector(n,init_q));




		// Build the way boxes will be bisected.
		// "LargestFirst" means that the dimension bisected
		// is always the largest one.
		LargestFirst lf;
		list<IntervalVector> s;
		s.push_back(init_box);
		IntervalVector box(3*n);
		cout<<"start"<<endl;
		int i=0;
		while (!s.empty()) {
			box=s.front();
			s.pop_front();
			contract_and_draw(outside,box,"[blue]k","[magenta]k");
			if (box.is_empty())  continue;

			contract_and_draw(inside,box,"[red]k","[cyan]k");
			if (box.is_empty())  continue;

			if (box.max_diam()<epsilon) {
				vibes::drawBox(box[ind].lb(),box[ind].ub(),box[ind+1].lb(),box[ind+1].ub(),"[yellow]k");
			} else {
				pair<IntervalVector,IntervalVector> boxes=lf.bisect(box);
				s.push_back(boxes.first);
				s.push_back(boxes.second);
			}
		}
	}

	vibes::endDrawing();      // <== closes the VIBES "connection"
	return 0;
}

