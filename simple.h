#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>

#include "cutility.h"
#include "MersenneTwister.h"

//#include <math.h>
using namespace std;

#ifndef Pi
#define Pi 3.141592653589793238462643
#endif

// Functions

double max(double num1, double num2);

int round (double x);

void VectorAdd(const int Flag, const vector<double>& x, const vector<double>& y, vector<double>& z);

void MMV(const vector<vector<double> >& A, const vector<double>& b, vector<double>& c);

double VectorNorm1(vector<double>& x);

double CND( double X );

double LinearInterp(vector<double>& xVec, vector<double>& yVec, double x);

double BlackScholes(char CallPutFlag, double S, double X, double T, double r, double v);

double EuroCallMC(double S, double K, double T, double r, double sigma, long N);

double EuroLattice(int Flag, double S, double K, double T, double r, double sigma, long N);

double AmLattice(int Flag, double S, double K, double T, double r, double sigma, long N);

double EuroCallExplicit(double S0, double K, double T, double r, double sigma, unsigned short M, unsigned short N);

double AmPutExplicit(double S0, double K, double T, double r, double sigma, unsigned short M, unsigned short N);

double AmPutCN(double S0, double K, double r, double T, double sigma, unsigned short M, unsigned short N, double omega, double tol);

double ImpVol(double F0, double K, double T, double r, double value, double v1, double v2, double tol, int N);

double f(char flag, double F0, double K, double T, double r, double sigma, double value);

// Implied Volatility (Newton Raphson)

double ImpVol2(double F0, double K, double T, double r, double value, double estimate, double tol, int N);

// Implied Volatlity ("proper" bisection)

double ImpVol3(double F0, double K, double T, double r, double value, double v1, double v2, double tol, int N);

// Return Discount

double ReturnDiscount(double * tVec, double * rVec, int N, double t);

// reverses an array of doubles

void reverse(double *A, int nA);

// Euler's constant

double Q29();