#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>

#include "cutility.h"
#include "MersenneTwister.h"

//#include <math.h>
using namespace std;

#ifndef Pi
#define Pi 3.141592653589793238462643
#endif

bool isNumeric(string stringToCheck); // string or a number

double linearInterp(const double& thisDay, const map<double, double>& crv); // linear interpolation

double timeWeightedRate(const double& thisDay, const map<double, double>& crv); // time weighted linear interpolation

double getForwardRate(const double& start, const double& end, const map<double, double>& crv); // forward rate

double getDF(const double& thisDay, const map<double, double>& crv); // time weighted linear rates, and the discount factor

double timeWeightedVariance(const double& thisDay, const map<double, double>& vol); // time weighted linear interploation of variances

double getForwardVar(const double& start, const double& end, const map<double, double>& vol); // forward variances

// Dot Product

double DotProduct(const vector<double>& x, const vector<double>& y);

// Matrix Transpose

void MatrixTranspose(const vector<vector<double> >& A, vector<vector<double> >& B);

// Matrix Multiplication

void MatrixMult(const vector<vector<double> >& A, const vector<vector<double> >& B, vector<vector<double> > & C);

// Time weighted interpolation

double TimeWeightedLinearInterp(const map<double, double>& data, const double t);

// Cholesky Decomposition

void Cholesky(const vector<vector<double> >& corr, vector<vector<double> >& A);

// Rainbow Option

double Rainbow(const vector<double>& weights, const vector<map<double, double> >& rates, const  map<double, double>& discounts,
			const vector<map<double, double> >& vols, const vector<double>& fixings,
			const vector<double>& spots, const vector<vector<double> >& corr, double T, int M);

// Asian Option

double Asian(int flag, double K, double T, const vector<map<double, double> >& rates, const vector<map<double, double> >& vols,
			 const map<double, double>& discounts, const vector<double>& weights, const vector<double>& fixings,
			const vector<double>& times, const vector<double>& spots, const vector<vector<double> >& corr, int N);

// Read Data (Rainbow)

bool GetData(vector<double>& weights, vector<map<double, double> >& rates, map<double, double>& discounts,
			 vector<map<double, double> >& vols, vector<double>& fixings,
			 vector<double>& spots, vector<vector<double> >& corr);