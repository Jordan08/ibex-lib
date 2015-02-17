//============================================================================
//                                  I B E X
// File        : exo04.cpp
// Author      : Jordan Ninin
// License     : See the LICENSE file
// Created     : Jan 29, 2014
// Last Update : Jan 29, 2014
//============================================================================


#include "ibex.h"

using namespace std;
using namespace ibex;



int main() {
	{
//		std::string figureName = "SIVIA_0";
//		vibes::beginDrawing();           // <== Initializes the VIBES "connection"
//		vibes::newFigure(figureName);       // <== Create a new VIBes figure

		double epsilon = 1.e-8;
		double gaol_prec= 1.e-6;
		double epsilon_time = 1.e-4;
		double T_final= 1;
		double time_out= 1000;

		Interval init_time(0,T_final);

		double secu =0.01;
		// number of plane
		int n=4;
		
		// initial position and speed
		Matrix P(n,2); // initial position of each plane
		Matrix V(n,2); // initial speed vector of each plane

		P[0][0] = 0;
		P[0][1] = 0;
		V[0][0] = 1;
		V[0][1] = 1;

		P[1][0] = 1;
		P[1][1] = 0;
		V[1][0] = -1;
		V[1][1] = 1;

		P[2][0] = 0;
		P[2][1] = 1;
		V[2][0] = 1.3;
		V[2][1] = -1;

		P[3][0] = 1;
		P[3][1] = 1;
		V[3][0] = -1;
		V[3][1] = -1.3;


		
		// Symbolic part
		Variable t,tm, delta,p(2),v(2),q, Z(2);
		Variable Tm(n),Delta(n),Q(n);

		Function traj_avion(t,tm,delta,q, p,v,
		  Return(
			chi(t-tm,
				p[0]+v[0]*t,
				chi(t-(tm+delta),
						p[0]+(1-q)*v[0]*tm+q*v[0]*t,
						p[0]+(q-1)*v[0]*delta + v[0]*t
				  )
			),
			chi(t-tm,
				p[1]+v[1]*t,
				chi(t-(tm+delta),
						p[1]+(1-q)*v[1]*tm+q*v[1]*t,
						p[1]+(q-1)*v[1]*delta + v[1]*t
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
		CtcForAll outside(outside1,mask,IntervalVector(1,init_time),epsilon_time);


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
		CtcExist inside(inside1,mask,IntervalVector(1,init_time),epsilon_time);

		// Build the initial box of V
		IntervalVector init_box(3*n,init_time);
		init_box.put(2*n,IntervalVector(n,Interval(-1.1,1.1)));


		// objective function
		const ExprNode* e=&(sqr(Q[0]));
		for (int i=1; i<n; i++)
			e = & (*e + sqr(Q[i]));
		Function norm_n(Q,*e);
		Function f_cost(Tm,Delta,Q, norm_n(Q)+norm_n(Delta));

		// create a stategy to bisect a element
		LargestFirst  bsc(epsilon);
		
		OptimCtc o(outside,inside,f_cost,bsc,epsilon,gaol_prec,gaol_prec);
		// the trace
			o.trace=1;

			// the allowed time for search
			o.timeout=time_out;

			// the search itself
			o.optimize(init_box);

			// printing the results
			o.report();
	}
//	 vibes::endDrawing();      // <== closes the VIBES "connection"
	return 0;
}

