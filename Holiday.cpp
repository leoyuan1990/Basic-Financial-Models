#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>

#include "cutility.h"
#include "MersenneTwister.h"
#include "simple.h"
#include "LessSimple.h"

//#include <math.h>
using namespace std;

#ifndef Pi
#define Pi 3.141592653589793238462643
#endif

int main (int argc, const char* argv[] )
{
	char key;
	int Flag;
	double S0;
	double K;
	double T;
	double r;
	double sigma;
	double price;
	double c;
	double F0;
	double impliedVol;
	int N;
	int i;
	int j;
		
	S0 = 50;
	K = 50;
	T = 5;
	r = 0.04;
	sigma = 0.25;
	N = 1000000;

	//price = EuroCallMC(S0,K,T,r,sigma,N);
	//cout << "EuroCallMC\t" << price << endl;

	//N = 30;
	//Flag = 1;
	//price = EuroLattice(Flag,S0,K,T,r,sigma,N);
	//cout << "EuroCallLattice\t" << price << endl;

	price = BlackScholes('c',S0,K,T,r,sigma);
	cout << "EuroCallBS\t" << price << endl;
	cin >> r;
	
	//price = AmLattice(Flag,S0,K,T,r,sigma,N);
	//cout << "AmCallLattice\t" << price << endl;

	//Flag = -1;
	//price = EuroLattice(Flag,S0,K,T,r,sigma,N);
	//cout << "EuroPutLattice\t" << price << endl;

	//price = AmLattice(Flag,S0,K,T,r,sigma,N);
	//cout << "AmPutLattice\t" << price << endl;

	//double M = 50;
	//N = 500;
	////price = EuroCallExplicit(S0,K,T,r,sigma,M,N);
	////cout << "EuroCallExp\t" << price << endl;

	//S0 = 50;
	//K = 50;
	//r = 0.1;
	//T = 5.0/12;
	//sigma = 0.4;

	//double omega = 1.2;
	//double tol = 1e-3;

	//price = AmPutCN(S0, K, r, T, sigma, M, N, omega, tol);
	//cout << "AmPutCN\t" << price << endl;
	//price = AmPutExplicit(S0, K, T, r, sigma, M, N);
	//cout << "AmPutExp\t" << price << endl;

	//F0 = 61;
	//T = 5;
	//r = 0.04;
	//K = 50;
	//c = 15;
	//double low = 0.1;
	//double high = 0.4;
	//double tolerance = 0.000001;
	//int lim = 100;	// max number of iterations

	//impliedVol = ImpVol(F0,K,T,r,c,low,high,tolerance,lim);
	//cout << "ImpliedVol1\t" << impliedVol << endl;
	//price = BlackScholes('c',F0*exp(-r*T),K,T,r,impliedVol);
	//cout << "ImpliedPrice1\t" << price << endl;

	//double guess = 0.3;

	//impliedVol = ImpVol2(F0,K,T,r,c,guess,tolerance,lim);
	//cout << "ImpliedVol2\t" << impliedVol << endl;
	//price = BlackScholes('c',F0*exp(-r*T),K,T,r,impliedVol);
	//cout << "ImpliedPrice2\t" << price << endl;

	//impliedVol = ImpVol3(F0,K,T,r,c,low,high,tolerance,lim);
	//cout << "ImpliedVol3\t" << impliedVol << endl;
	//price = BlackScholes('c',F0*exp(-r*T),K,T,r,impliedVol);
	//cout << "ImpliedPrice3\t" << price << endl;

	//double discfact;
	//double timeVec[] = {0.25, 0.50, 1.00, 1.50};
	//double rateVec[] = {0.02, 0.03, 0.04, 0.05};
	//double time = 270.0/365;
	//int arraySize = 4;
	//discfact = ReturnDiscount(timeVec, rateVec, arraySize, time);
	//cout << "ReturnDiscount\t" << discfact << endl; 
	
	//map<double, double> data;
	//data.insert(make_pair(0.25, 0.02));
	//data.insert(make_pair(0.50, 0.03));
	//data.insert(make_pair(1.00, 0.04));
	//data.insert(make_pair(1.50, 0.05));

	//double rate = TimeWeightedLinearInterp(data,time);
	//cout << "TimeWeightedLinearInterp\t" << rate << endl;
	
	//vector<vector<double> > A;
	//vector<double> x;
	//vector<double> b;
	//
	//// resizes the vectors and matrices

	//x.resize(M-1);
	//b.resize(M-1);
	//A.resize(M-1);
	//for (i = 0; i < M-1; ++i)
	//{
	//	A[i].resize(M-1);
	//}

	//vector<double> weights;
	//vector<map<double, double> > rates;
	//map<double, double> discounts;
	//vector<map<double, double> > vols;
	//vector<double> fixings, spots;
	//vector<vector<double> > corr;
	//vector<vector<double> > A, B, C;

	//GetData(weights, rates, discounts, vols, fixings, spots, corr);

	////for(i=0; i<corr.size(); ++i)
	////{
	////	for(j=0; j<corr[i].size(); ++j)		// NEVER use j<corr[0], always use j<corr[i]
	////		cout << corr[i][j] << "\t";

	////	cout << endl;
	////}

	////Cholesky(corr, A);
	////MatrixTranspose(A, B);
	////MatrixMult(A, B, C);

	////for(i=0; i<C.size(); ++i)
	////{
	////	for(j=0; j<C[i].size(); ++j)	
	////		cout << C[i][j] << "\t";

	////	cout << endl;
	////}

	//T = 6;
	//N = 100000;
	//Flag = 1;
	//vector<double> times;
	//times.push_back(1.5);
	//times.push_back(2.5);
	//times.push_back(3.5);
	//times.push_back(4.5);
	//times.push_back(5.5);
	//
	///*price = Rainbow(weights, rates, discounts, vols, fixings, spots, corr, T, N);
	//cout << "Rainbow\t" << price << endl;*/

	//K = 0.8;

	//price = Asian(Flag, K, T, rates, vols, discounts, weights, fixings,	times, spots, corr, N);
	//cout << "Asian\t" << price << endl;
	//
	//cin >> key;

	/*N = 6;
	double B [6] = {3, 4, 2, 8, 19, 20};
	reverse(B, N);
	for (i=0; i<N; ++i)
	{
		cout << B[i] << endl;
	}*/

	//int n = Q29();
	//cout << n << endl;

	return 0;
}