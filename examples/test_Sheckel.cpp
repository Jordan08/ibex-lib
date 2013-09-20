//============================================================================
//                                  I B E X
// File        : arith03.cpp
// Author      : Jordan Ninin
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Fev 28, 2013
// Last Update : Fev 28, 2013
//============================================================================

#include "ibex.h"
#include <time.h>

using namespace std;
using namespace ibex;

int main() {


	{
		cout << "==========================================" << endl;
		cout << "==========================================" << endl;
		int n = 100000;
		cout << "TEST 3 : " << n << " evaluations of the Sheckel-5 Function "<< endl;
		double A[5][4] = { { 4, 4, 4, 4 }, { 1, 1, 1, 1 }, { 8, 8, 8, 8 }, { 6,
				6, 6, 6 }, { 3, 7, 3, 7 } };

		double c[5] = { 0.1, 0.2, 0.2, 0.4, 0.4 };
		clock_t start, endtime;
		double cpuTime;
		{
			double f, z;
			double x[4] = { 4, 4, 4, 4 };
			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					for (int j = 0; j < 4; j++) {
						z = z + pow((x[j] - A[i][j]), 2);
					}
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / ((double) CLOCKS_PER_SEC);
			cout << "double : CPU-time = " << cpuTime << " seconds" << endl;
		}
		{
			Interval f, z;
			IntervalVector x(4, Interval(3.9, 4.1));
			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					for (int j = 0; j < 4; j++) {
						z = z + pow((x[j] - A[i][j]), 2);
					}
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Interval : CPU-time = " << cpuTime << " seconds" << endl;
		}
		{

			// Default is AF_fAF2
			Affine2 f, z;
			Affine2Vector x(4, Interval(3.9, 4.1), true );  // Initialization with x[i] = Affine2(4,i+1,Interval(3.9, 4.1));
			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					for (int j = 0; j < 4; j++) {
						z = z + pow((x[j] - A[i][j]), 2);
					}
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine2 Default : CPU-time = " << cpuTime << " seconds" << endl;
		}

		{
			Affine2Main<AF_fAF1> f, z;
			Interval initIA (3.9, 4.1);
			Affine2Main<AF_fAF1> x1 (4, 1, initIA);
			Affine2Main<AF_fAF1> x2 (4, 1, initIA);
			Affine2Main<AF_fAF1> x3 (4, 1, initIA);
			Affine2Main<AF_fAF1> x4 (4, 1, initIA);

			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					z = z + pow((x1 - A[i][0]), 2);
					z = z + pow((x2 - A[i][1]), 2);
					z = z + pow((x3 - A[i][2]), 2);
					z = z + pow((x4 - A[i][3]), 2);
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine2 AF_fAF1 : CPU-time = " << cpuTime << " seconds" << endl;
		}


		{
			Affine2Main<AF_fAF2> f, z;
			Interval initIA (3.9, 4.1);
			Affine2Main<AF_fAF2> x1 (4, 1, initIA);
			Affine2Main<AF_fAF2> x2 (4, 1, initIA);
			Affine2Main<AF_fAF2> x3 (4, 1, initIA);
			Affine2Main<AF_fAF2> x4 (4, 1, initIA);

			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					z = z + pow((x1 - A[i][0]), 2);
					z = z + pow((x2 - A[i][1]), 2);
					z = z + pow((x3 - A[i][2]), 2);
					z = z + pow((x4 - A[i][3]), 2);
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine2 AF_fAF2 : CPU-time = " << cpuTime << " seconds" << endl;
		}


		{
			Affine2Main<AF_fAF2_fma> f, z;
			Interval initIA (3.9, 4.1);
			Affine2Main<AF_fAF2_fma> x1 (4, 1, initIA);
			Affine2Main<AF_fAF2_fma> x2 (4, 1, initIA);
			Affine2Main<AF_fAF2_fma> x3 (4, 1, initIA);
			Affine2Main<AF_fAF2_fma> x4 (4, 1, initIA);

			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					z = z + pow((x1 - A[i][0]), 2);
					z = z + pow((x2 - A[i][1]), 2);
					z = z + pow((x3 - A[i][2]), 2);
					z = z + pow((x4 - A[i][3]), 2);
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine2 AF_fAF2_fma : CPU-time = " << cpuTime << " seconds" << endl;
		}

		{
			Affine2Main<AF_iAF> f, z;
			Interval initIA (3.9, 4.1);
			Affine2Main<AF_iAF> x1 (4, 1, initIA);
			Affine2Main<AF_iAF> x2 (4, 1, initIA);
			Affine2Main<AF_iAF> x3 (4, 1, initIA);
			Affine2Main<AF_iAF> x4 (4, 1, initIA);

			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					z = z + pow((x1 - A[i][0]), 2);
					z = z + pow((x2 - A[i][1]), 2);
					z = z + pow((x3 - A[i][2]), 2);
					z = z + pow((x4 - A[i][3]), 2);
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine2 AF_iAF : CPU-time = " << cpuTime << " seconds" << endl;
		}
		{
			Affine2Main<AF_No> f, z;
			Interval initIA (3.9, 4.1);
			Affine2Main<AF_No> x1 (4, 1, initIA);
			Affine2Main<AF_No> x2 (4, 1, initIA);
			Affine2Main<AF_No> x3 (4, 1, initIA);
			Affine2Main<AF_No> x4 (4, 1, initIA);

			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					z = z + pow((x1 - A[i][0]), 2);
					z = z + pow((x2 - A[i][1]), 2);
					z = z + pow((x3 - A[i][2]), 2);
					z = z + pow((x4 - A[i][3]), 2);
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine2 AF_No : CPU-time = " << cpuTime << " seconds" << endl;
		}

		/*
		 * ==========================================
		 * with -O3 option
		 * ./waf configure --prefix=/home/nininjo/Documents/WORK/PROJET/IBEX/OUT --with-soplex=/home/nininjo/Logiciel/soplex-1.7.1 --with-filib=/home/nininjo/Logiciel/filib
		 * ./waf
		 * __build__/examples/test_Sheckel
TEST 3 : 100000 evaluations of the Sheckel-5 Function
double : CPU-time = 0 seconds
Interval : CPU-time = 0.14 seconds
Affine2 Default : CPU-time = 1.49 seconds
Affine2 AF_fAF1 : CPU-time = 1.79 seconds
Affine2 AF_fAF2 : CPU-time = 1.43 seconds
Affine2 AF_fAF2_fma : CPU-time = 2.68 seconds
Affine2 AF_iAF : CPU-time = 3.17 seconds
Affine2 AF_No : CPU-time = 1.31 seconds
		 *
		 * ==========================================
		 * without optimization
		 * ./waf configure --prefix=/home/nininjo/Documents/WORK/PROJET/IBEX/OUT --with-soplex=/home/nininjo/Logiciel/soplex-1.7.1 --with-filib=/home/nininjo/Logiciel/filib --with-debug
		 * ./waf
		 * __build__/examples/test_Sheckel
TEST 3 : 100000 evaluations of the Sheckel-5 Function
double : CPU-time = 0.03 seconds
Interval : CPU-time = 2.26 seconds
Affine2 Default : CPU-time = 22.61 seconds
Affine2 AF_fAF1 : CPU-time = 29.18 seconds
Affine2 AF_fAF2 : CPU-time = 22.16 seconds
Affine2 AF_fAF2_fma : CPU-time = 22.72 seconds
Affine2 AF_iAF : CPU-time = 63.33 seconds
Affine2 AF_No : CPU-time = 19.19 seconds
		 *
		 *
		 *
		 */


	}

	return 0;
}

