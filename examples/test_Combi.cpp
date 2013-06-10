#include "ibex.h"

using namespace std;
using namespace ibex;



int main(int argc, char** argv) {

	try {

		// Load a system of equations
		System sys(argv[1]);

		// the extended system
		System ext_sys(sys, System::EXTEND);

		cout << " " << argv[1] << " &" ;

		double prec       = 1.e-6;
		double goal_prec  = 1.e-6;
		double time_limit = 600;


		if (!sys.goal) {
			ibex_error(
					" input file has not goal (it is not an optimization problem).");
		}

		Bsc * bs = new SmearMaxRelative(ext_sys, prec);

		// The contractors

		// the first contractor called
		CtcHC4 hc4(ext_sys.ctrs, 0.01, true);
		// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
		CtcHC4 hc44cid(ext_sys.ctrs, 0.1, true);
		// hc4 inside xnewton loop
		CtcHC4 hc44xn(ext_sys.ctrs, 0.01, false);

		// The ACID contractor (component of the contractor  when filtering == "acidhc4")
		CtcAcid acidhc4(ext_sys, hc44cid, true);
		// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4"
		CtcCompo hc4acidhc4(hc4, acidhc4);


		// the linear relaxation contractor: ART
		CtcARTiter ctcart(ext_sys);
		// ART linear relaxation , hc4  with default fix point ratio 0.2
		CtcART ca(ctcart, hc44xn);
		//  the actual contractor  ctc + ART
		CtcCompo ctc_art(hc4acidhc4, ca);
		// one point probed when looking for a new feasible point (updating the loup)


		// The CtcXNewtonIter contractor
		// corner selection for linearizations : two corners are slected, a random one and its opposite
		vector<CtcXNewtonIter::corner_point> cpoints;
		cpoints.push_back(CtcXNewtonIter::RANDOM);
		cpoints.push_back(CtcXNewtonIter::RANDOM_INV);

		// the linear relaxation contractor
		CtcXNewtonIter ctcxnewton (ext_sys,cpoints);

		// fixpoint linear relaxation , hc4  with default fix point ratio 0.2
		CtcXNewton cxn (ctcxnewton, hc44xn);
		//  the actual contractor  ctc + xnewton
		CtcCompo ctc_xn (hc4acidhc4, cxn);

		CtcCompo ctc_compo(ctc_xn,ctc_art);

		int samplesize = 1;
		// the optimizer : the same precision goal_prec is used as relative and absolute precision
		Optimizer o(sys, *bs, ctc_compo, prec, goal_prec, goal_prec, samplesize);

		// This option limits the search time
		o.timeout=time_limit;

		// This option prints each better feasible point when it is found
		o.trace=0;

		// display solutions with up to 12 decimals
		cout.precision(12);

		// Search for the optimum
		o.optimize(sys.box);

		// Report some information (computation time, etc.)
		//o.report();
		o.report_perf();


		return 0;

	}
	catch(ibex::SyntaxError& e) {
		cout << e << endl;
	}
}
