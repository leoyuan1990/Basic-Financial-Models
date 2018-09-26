#include "simple.h"

// function returning the max between two numbers
 
double max(double num1, double num2) 

{
   // local variable declaration
   double result;
 
   if (num1 > num2)
      result = num1;
   else
      result = num2;
 
   return result; 
}

// Round function

int round (double x)

{
	if (x - floor(x) < 0.5)
		return floor(x);
	else
		return floor(x)+1;
}

// Vector Addition

void VectorAdd(const int Flag, const vector<double>& x, const vector<double>& y, vector<double>& z)
{
	
	double xSize = x.size();
	double ySize = y.size();
	z.clear();
	z.resize(xSize);

	if (xSize != ySize)
	{
		cout << "Vector sizes must be equal." << endl;
		exit(EXIT_FAILURE);
	}
	else if (Flag != -1 && Flag != 1)
	{
		cout << "Flag must equal -1 for subtraction or 1 for addition." << endl;
		exit(EXIT_FAILURE);
	}

	for(int i=0; i<xSize; ++i)
	{
		z[i] = x[i] + Flag * y[i];
	}

}

// Matrix Multiply Vector

void MMV(const vector<vector<double> >& A, const vector<double>& b, vector<double>& c)	// need a space after double (ie > > NOT >>)
{

	double Asize = A.size();
	double bSize = b.size();
	c.resize(Asize);

	for(int i=0; i<A.size(); ++i)
	{
		c[i] = 0.0;
		for(int j=0; j<b.size(); ++j)
			c[i] += A[i][j] * b[j];
	}
		
}

// Vector Norm (L1)

double VectorNorm1(vector<double>& x)
{
	double answer = 0;
	double size = x.size();

	for (int i=0; i<size; i++)
	{
		answer = answer + fabs(x[i]);
	}
	
	return answer;

}

// The cumulative normal distribution function
double CND( double X )

{
	double L, K, w ;

	double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
	double const a4 = -1.821255978, a5 = 1.330274429;

	L = fabs(X);
	K = 1.0 / (1.0 + 0.2316419 * L);
	w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

	if (X < 0 ){
		w= 1.0 - w;
	}
	return w;
}

// linear interpolation

double LinearInterp(vector<double>& xVec, vector<double>& yVec, double x)
{

	int xsize = xVec.size();
	int ysize = yVec.size();
	if (xsize != ysize)
	{
		cout << "The vector sizes must be equal" << endl;
		exit(EXIT_FAILURE);
	}
	int i=0;
	int N = xsize;
	double x1;
	double x2;
	double y1;
	double y2;
	double slope;

	
	if (x >= xVec[0] && x <= xVec[N-1])
	{
		while (xVec[i] < x)
		{
			++i;
		}
		x1 = xVec[i-1];
		x2 = xVec[i];
		y1 = yVec[i-1];
		y2 = yVec[i];
		slope = (y2-y1)/(x2-x1);
		return y1+slope*(x-x1);
	}
	else if (x < xVec[0])
	{
		slope = (yVec[1]-yVec[0])/(xVec[1]-xVec[0]);
		return yVec[0]+slope*(x-xVec[0]);
	}
	else
	{
		slope = (yVec[N-1]-yVec[N-2])/(xVec[N-1]-xVec[N-2]);
		return yVec[N-1]+slope*(x-xVec[N-1]);
	}
}

// The Black and Scholes (1973) Stock option formula
double BlackScholes(char CallPutFlag, double S, double X, double T, double r, double v)
{
	double d1, d2;

	d1=(log(S/X)+(r+v*v/2)*T)/(v*sqrt(T));
	d2=d1-v*sqrt(T);

	if(CallPutFlag == 'c')
		return S *cutility::svCumNorm(d1)-X * exp(-r*T)*cutility::svCumNorm(d2);
	else if(CallPutFlag == 'p')
		return X * exp(-r * T) * cutility::svCumNorm(-d2) - S * cutility::svCumNorm(-d1);
	else
		exit(EXIT_FAILURE);

}

// Monte Carlo European Call pricing function

double EuroCallMC(double S, double K, double T, double r, double sigma, long N)

{
	double U1, U2;
	double X;
	double nuT;
	double siT;
	double price = 0.0;	// how to change this into an array?
	double discount = exp(-r*T);
	long i = 0;

	init_genrand(static_cast<unsigned long>(time(NULL))); 

	nuT = (r-0.5*sigma*sigma)*T; 
	siT = sigma*sqrt(T);

	for(i=0; i<N; ++i){
		X = gasdev();
		price += max(0, S*exp(nuT+siT*X)-K)*discount;
	}

	return price/N;
}

// Lattice Method (European)

double EuroLattice(int Flag, double S, double K, double T, double r, double sigma, long N)

{
	if (Flag != 1 && Flag != -1)
	{
		cout << "Enter 1 or -1" << endl;
		return -9999;
	}

	double deltaT = T/N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1/u;
	double p = (exp(r*deltaT)-d)/(u-d);
	double q = 1-p;
	vector<double> S_Val(2*N+1);
	vector<double> C_Val(2*N+1);
	double discount;
	long i;
	long tau;

	discount = exp(-r*deltaT);
	S_Val[0] = S*pow(d,N);
	double p_u = p*discount;
	double p_d = q*discount;
	for (i=2;i<2*N+1;i+=2)
	{
		S_Val[i] = u*u*S_Val[i-2];
	}

	for (i=0;i<2*N+1;i+=2)
	{
		C_Val[i] = max(Flag*(S_Val[i]-K),0);
	}

	for (tau=0;tau<N;++tau)
	{
		for (i=tau+1;i<2*N-tau;i+=2)
		{
			C_Val[i] = p_u*C_Val[i+1]+p_d*C_Val[i-1];
		}
	}
		
	return C_Val[N];
}

// Lattice Method (American)

double AmLattice(int Flag, double S, double K, double T, double r, double sigma, long N)

{
	if (Flag != 1 && Flag != -1)
	{
		cout << "Enter 1 or -1" << endl;
		return -9999;
	}
	
	double deltaT = T/N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1/u;
	double p = (exp(r*deltaT)-d)/(u-d);
	double q = 1-p;
	vector<double> S_Val(2*N+1);
	vector<double> C_Val(2*N+1);
	double discount;
	long i;
	long tau;
	/*sort(S_Val.begin(),S_Val.end());*/
	discount = exp(-r*deltaT);
	S_Val[0] = S*pow(d,N);
	double p_u = p*discount;
	double p_d = q*discount;
	for (i=1;i<2*N+1;++i)
	{
		S_Val[i] = u*S_Val[i-1];
	}

	for (i=0;i<2*N+1;++i)
	{
		C_Val[i] = max(Flag*(S_Val[i]-K),0);
	}

	double hold;
	double exercise;

	for (tau=0;tau<N;++tau)
	{
		for (i=tau+1;i<2*N-tau;i+=2)
		{
			hold = p_u*C_Val[i+1]+p_d*C_Val[i-1];
			exercise = Flag*(S_Val[i]-K);
			C_Val[i] = max(hold,exercise);
		}
	}
		
	return C_Val[N];
}

// Explicit FD

double EuroCallExplicit(double S0, double K, double T, double r, double sigma, unsigned short M, unsigned short N)

{
	double factor = 2;
	double Smax = factor*S0;
	double dS = Smax/M;
	double dt = T/N;
	int i;
	int j;

	//if (dt > 0.5 * dS * dS)
	//{
	//	cout << "Will not be stable" << endl;
	//	exit(EXIT_FAILURE);
	//}

	// solution matrix
	vector<vector<double>> sol;
	sol.resize(N+1);
	for(i=0; i<N+1; ++i)
	{
		sol[i].resize(M+1);
	}

	vector<double> vetS(M+1);
	vector<double> veti(M+1);
	vector<double> vetj(N+1);

	// coefficient vectors

	vector<double> a(M+1);
	vector<double> b(M+1);
	vector<double> c(M+1);

	// initializes the vectors

	for (i=0; i<M+1; ++i)
	{
		vetS[i] = i*dS;
		veti[i] = i;
		a[i] = 0.5*dt*(sigma*sigma*veti[i]-r)*veti[i];
		b[i] = 1-dt*(sigma*sigma*veti[i]*veti[i]+r);
		c[i] = 0.5*dt*(sigma*sigma*veti[i]+r)*veti[i];
	}

	for (j=0; j<N+1; ++j)
	{
		vetj[j]= j;
	}

	// boundary conditions

	for (j=0; j<N+1; ++j)
	{
		sol[j][0] = 0;
		sol[j][M] = Smax-K*exp(-r*(N-vetj[j])*dt);
	}

	// terminal condition

	for (i=0; i<M+1; ++i)
	{
		sol[N][i] = max(vetS[i]-K,0);
	}

	// solution

	for (j=N-1; j>-1; --j)
	{
		for (i=1; i<M; ++i)
		{
			sol[j][i] = a[i]*sol[j+1][i-1]+b[i]*sol[j+1][i]+c[i]*sol[j+1][i+1];
		}
	}

	return LinearInterp(vetS,sol[0],S0);
}

// Explicit FD Am Put

double AmPutExplicit(double S0, double K, double T, double r, double sigma, unsigned short M, unsigned short N)

{
	double factor = 2;
	double Smax = factor*S0;
	double dS = Smax/M;
	double dt = T/N;
	int i;
	int j;

	//if (dt > 0.5 * dS * dS)
	//{
	//	cout << "Will not be stable" << endl;
	//	exit(EXIT_FAILURE);
	//}

	// solution matrix
	vector<vector<double>> sol;
	sol.resize(N+1);
	for(i=0; i<N+1; ++i)
	{
		sol[i].resize(M+1);
	}

	vector<double> vetS(M+1);
	vector<double> veti(M+1);
	vector<double> vetj(N+1);

	// coefficient vectors

	vector<double> a(M+1);
	vector<double> b(M+1);
	vector<double> c(M+1);

	// initializes the vectors

	for (i=0; i<M+1; ++i)
	{
		vetS[i] = i*dS;
		veti[i] = i;
		a[i] = 0.5*dt*(sigma*sigma*veti[i]-r)*veti[i];
		b[i] = 1-dt*(sigma*sigma*veti[i]*veti[i]+r);
		c[i] = 0.5*dt*(sigma*sigma*veti[i]+r)*veti[i];
	}

	for (j=0; j<N+1; ++j)
	{
		vetj[j]= j;
	}

	// boundary conditions

	for (j=0; j<N+1; ++j)
	{
		sol[j][0] = K*exp(-r*(N-vetj[j])*dt);
		sol[j][M] = 0;
	}

	// terminal condition

	for (i=0; i<M+1; ++i)
	{
		sol[N][i] = max(K-vetS[i],0);
	}

	// solution

	for (j=N-1; j>-1; --j)
	{
		for (i=1; i<M; ++i)
		{
			sol[j][i] = a[i]*sol[j+1][i-1]+b[i]*sol[j+1][i]+c[i]*sol[j+1][i+1];
			sol[j][i] = max(sol[j][i], K - i * dS);
		}
	}

	return LinearInterp(vetS,sol[0],S0);
}

// AmPutCN

double AmPutCN(double S0, double K, double r, double T, double sigma, unsigned short M, unsigned short N, double omega, double tol)
{

	double factor = 2;
	double Smax = factor*S0;
	double dS = Smax/M;
	double dt = T/N;
	double realmax = 1e+15;
	double error;
	int i;
	int j;
	
	vector<double> oldVal(M-1);
	vector<double> newVal(M-1);
	vector<double> pastVal(M-1);
	vector<double> differenceVal(M-1);
	vector<double> boundVal(N+1);
	vector<double> payoff(M-1);
	
	// rhs
	double aux;
	vector<double> rhs(M-1);
	vector<double> M2timesF;

	vector<double> vetS(M+1);
	vector<double> veti(M+1);
	vector<double> vetj(N+1);
	
	vector<vector<double> > M2;
	M2.resize(M-1);
	for (i=0; i<M-1; i++)
	{
		M2[i].resize(M-1);
	}

	// coefficient vectors
	vector<double> a(M+1);
	vector<double> b(M+1);
	vector<double> c(M+1);

	vetS[0] = 0;
	veti[0] = 1;
	a[0] = 0.25*dt*(sigma*sigma*veti[0]-r)*veti[0];
	b[0] = -0.5*dt*(sigma*sigma*veti[0]*veti[0]+r);
	c[0] = 0.25*dt*(sigma*sigma*veti[0]+r)*veti[0];

	for (i=1; i<M+1; ++i)
	{
		vetS[i] = vetS[i-1] + dS;
		veti[i] = veti[i-1] + 1;
		a[i] = 0.25*dt*(sigma*sigma*veti[i]-r)*veti[i];
		b[i] = -0.5*dt*(sigma*sigma*veti[i]*veti[i]+r);
		c[i] = 0.25*dt*(sigma*sigma*veti[i]+r)*veti[i];
	}

	vetj[0] = 0;
	for (j=1; j<N+1; ++j)
	{
		vetj[j]= vetj[j-1] + 1;
	}

	// terminal and boundary conditions

	for (i=0; i<M-1; ++i)
	{
		payoff[i] = max(K-vetS[i+1],0);
		pastVal[i] = payoff[i];
	}

	for (j=0; j<N+1; ++j)
	{
		boundVal[j] = K*exp(-r*dt*(N-vetj[j]));
	}

	// right hand side set up

	for (i=0; i<M-1; ++i)
	{
		for (j=0; j<M-1; ++j)
		{
			if (i == j + 1)
				M2[i][j] = a[i];
			else if (i == j)
				M2[i][j] = 1 + b[i];
			else if (i == j - 1)
				M2[i][j] = c[i];
		}
	}

	for (j=N-1; j>-1; --j)
	{
		aux = a[0] * (boundVal[j] + boundVal[j+1]);
		MMV(M2, pastVal, M2timesF);
		rhs = M2timesF;
		rhs[0] += aux;
		oldVal = pastVal;
		error = realmax;

		// iteration

		while (tol < error)
		{
			newVal[0] = max(payoff[0], oldVal[0] + omega/(1 - b[0]) * (rhs[0] - (1 - b[0]) * oldVal[0] 
			+ c[0] * oldVal[1]));

			for (int k=1; k<M-2; ++k)
			{
				newVal[k] = max(payoff[k], oldVal[k] + omega/(1 - b[k]) * (rhs[k] + a[k] * newVal[k-1] 
				- (1 - b[k]) * oldVal[k] + c[k] * oldVal[k+1]));
			}

			newVal[M-2] = max(payoff[M-2], oldVal[M-2] + omega/(1 - b[M-2]) * (rhs[M-2] + a[M-2] * newVal[M-3] 
			- (1 - b[M-2]) * oldVal[M-2]));
			
			VectorAdd(-1, newVal, oldVal, differenceVal);
			error = VectorNorm1(differenceVal);

			oldVal = newVal;
		}

		pastVal = newVal;

	}

	vector<double> answer(M+1);
	answer[0] = boundVal[0];
	answer[M] = 0;
	for (i=1; i<M; ++i)
	{
		answer[i] = newVal[i-1];
	}

	double final;
	final = LinearInterp(vetS, answer, S0);

	return final;
}

// Implied Volatility (using bisection method)

double ImpVol(double F0, double K, double T, double r, double value, double v1, double v2, double tol, int N)
{
	vector<double> lower(N);	// lower and upper BS prices using lower and upper estimates
	vector<double> upper(N);
	vector<double> sigma1(N);		// lower and upper estimates for volatility
	vector<double> sigma2(N);
	char flag = 'c';
	double sigmaTemp;
	int iterations = 0;

	lower[0] = BlackScholes(flag,F0*exp(-r*T),K,T,r,v1);	// initial lower volatility estimate
	upper[0] = BlackScholes(flag,F0*exp(-r*T),K,T,r,v2);	// initial upper volatility estimate

	// checks to see if appropriate initial estimates are entered
	if (v1 >= v2)
	{
		cout << "The lower estimate for the volatility must be lower than the upper estimate." << endl;
		exit(EXIT_FAILURE);
	}
	else if (lower[0] > value && upper[0] > value)
	{
		cout << "Lower volatility estimate is too high." << endl;
		exit(EXIT_FAILURE);
	}
	else if (upper[0] < value && lower[0] < value)
	{
		cout << "Upper volalitiy estimate is too low." << endl;
		exit(EXIT_FAILURE);
	}
	else if (lower [0] > value && upper[0] < value)
	{
		cout << "Lower volatility estimate is too high and upper volatility estimate is too low." << endl;
		exit(EXIT_FAILURE);
	}

	sigma1[0] = v1;
	sigma2[0] = v2;
	do
	{
		lower[iterations] = BlackScholes(flag,F0*exp(-r*T),K,T,r,sigma1[iterations]);
		upper[iterations] = BlackScholes(flag,F0*exp(-r*T),K,T,r,sigma2[iterations]);
		sigmaTemp = (sigma1[iterations] + sigma2[iterations])/2;
		if (fabs(value - lower[iterations]) < fabs(value - upper[iterations]))
		{
			sigma1[iterations + 1] = sigma1[iterations];
			sigma2[iterations + 1] = sigmaTemp;
		}
		else
		{
			sigma1[iterations + 1] = sigmaTemp;
			sigma2[iterations + 1] = sigma2[iterations];
		}
		iterations = iterations + 1;
	} while (iterations < N && (sigma2[iterations]-sigma1[iterations] > tol*(v2-v1)));

	return sigma2[iterations];
}

// Function (Newton Raphson)

double f(char flag, double F0, double K, double T, double r, double sigma, double value)
{
	return BlackScholes(flag,F0*exp(-r*T),K,T,r,sigma) - value;
}

// Implied Volatility (Newton Raphson)

double ImpVol2(double F0, double K, double T, double r, double value, double estimate, double tol, int N)
{

	double stepSize = 1e-8;		// step size used in forward approximation of derivative (h)
	double f1;
	double f2;				// f(x+h)
	double deriv;			// (f(x+h)-f(x))/h
	char flag = 'c';
	int iterations = 0;

	f1 = f(flag,F0*exp(-r*T),K,T,r,estimate,value);
	
	while (fabs(f1) > tol && iterations < N)
	{
		f1 = f(flag,F0,K,T,r,estimate,value);
		f2 = f(flag,F0,K,T,r,estimate + stepSize,value);
		deriv = (f2 - f1)/stepSize;
		estimate -= f1/deriv;
		iterations += 1;
	}
	
	return estimate;
}

// Implied Volatlity ("proper" bisection)

double ImpVol3(double F0, double K, double T, double r, double value, double v1, double v2, double tol, int N)
{

	double low = v1;
	double high = v2;
	double mid = (low+high)/2;
	int iterations = 0;
	char c = 'c';
	double lowF = f(c,F0,K,T,r,low,value);
	double midF = f(c,F0,K,T,r,mid,value);
	double highF = f(c,F0,K,T,r,high,value);

	if (lowF * highF > 0)
	{
		cout << "Initial estimates are inappropriate." << endl;
		exit(EXIT_FAILURE);
	}
	else if (lowF == 0)
		return low;
	else if (highF == 0)
		return high;

	while (iterations < N && fabs(midF) > tol)
	{
		if (lowF * midF < 0)
			high = mid;
		else
			low = mid;
		mid = (low + high)/2;
		lowF = f(c,F0,K,T,r,low,value);
		midF = f(c,F0,K,T,r,mid,value);
		highF = f(c,F0,K,T,r,high,value);
	}

	return mid;
}

// Return Discount

double ReturnDiscount(double * tVec, double * rVec, int N, double t)
{

	/*int tsize = tVec.size();
	int rsize = rVec.size();
	if (tsize != rsize)
	{
		cout << "The array sizes must be equal" << endl;
		exit(EXIT_FAILURE);
	}*/
	int i=0;
	double t1;
	double t2;
	double r1;
	double r2;
	double slope;
	double timeSlope;
	double discount;
	double final;
	
	if (t >= tVec[0] && t <= tVec[N-1])
	{
		while (tVec[i] < t)
		{
			++i;
		}
		t1 = tVec[i-1];
		t2 = tVec[i];
		r1 = rVec[i-1];
		r2 = rVec[i];
		//slope = (r2 - r1)/(t2 - t1);
		timeSlope = (t2 * r2 - t1 * r1)/(t2 - t1);
		final = t1*r1+timeSlope*(t - t1);
	}
	else if (t < tVec[0])
	{
		slope = (rVec[1] - rVec[0])/(tVec[1] - tVec[0]);
		timeSlope = (tVec[1] * rVec[1] - tVec[0] * rVec[0])/(tVec[1] - tVec[0]);
		final = tVec[0] * rVec[0] + timeSlope * (t - tVec[0]);
	}
	else
	{
		slope = (rVec[N-1] - rVec[N-2])/(tVec[N-1] - tVec[N-2]);
		timeSlope = (tVec[N-1] * rVec[N-1] - tVec[N-2] * rVec[N-2])/(tVec[N-1] - tVec[N-2]);
		final = tVec[N-1] * rVec[N-1] + timeSlope * (t - tVec[N-1]);
	}

	return (final/t);
	//return exp(-final);
}

// reverses an array of doubles

void reverse(double *A, int nA)

{
	int i;
	double tmp;

	if (nA % 2 == 0)
	{
		for (i=0; i<(nA - 2)/2; ++i)
		{
			tmp = A[i];
			A[i] = A[nA - i];
			A[nA] = tmp;
		}
	}
	else
		for (i=0; i<(nA - 3)/2; ++i)
		{
			tmp = A[i];
			A[i] = A[nA - i];
			A[nA] = tmp;
		}
	return;
}

// Euler's constant

double Q29()

{
	const int lim = 100000;
	double total = 0;
	double n = 2;

	while (fabs(total - 0.5772) >= 0.0001 && n < lim)
	{
		total = total + 1/n - log(n) + log(n - 1);
		++n;
	}

	return n;
}