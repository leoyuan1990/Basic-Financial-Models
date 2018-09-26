// gutility.cpp

#include "cutility.h"
#include <cmath>
#include <cstdlib>
#include <iostream>


bool cutility::leapYear(const long yyyy) 
{
	//A year is a leap year if it is divisible by 4, unless it is divisible by 100 
	//and not by 400. This means that 1984 and 2000 were leap years, but 1900 was not. 

	return (yyyy % 4 == 0 && yyyy % 100 != 0 || yyyy % 400 == 0);
}

int cutility::dayOfWeek(const long yyyymmdd)// 0 = sunday
{
   // From C FAQ, by Tomohiko Sakamoto
  int yyyy = yyyymmdd / 10000;
  int mm = (yyyymmdd - yyyy * 10000) / 100;
  int dd = yyyymmdd % 100;
  const int t[] = {0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4};
  yyyy -= mm < 3;

  return (yyyy + yyyy/4 - yyyy/100 + yyyy/400 + t[mm-1] + dd) % 7;
}


bool cutility::isWeekDay(const long yyyymmdd)
{
  int d = dayOfWeek(yyyymmdd);
  return (d>0 && d<6);
}

bool cutility::isWeekend(const long yyyymmdd)
{
  int d = dayOfWeek(yyyymmdd); // d=0: Sunday, d=6: Saturday.
  return (d==0 || d == 6);
}


long cutility::days30360( const long dt1, const long dt2)// dt is in format 'yyyymmdd'
{
	long d1 = gmin(dt1, dt2);
	long d2 = gmax(dt1, dt2);

	long yyyy1 = d1 / 10000;
	long mm1 = (d1 - yyyy1 * 10000) / 100;
	long dd1 = d1 % 100;

	long yyyy2 = d2 / 10000;
	long mm2 = (d2 - yyyy2 * 10000) / 100;
	long dd2 = d2 % 100;

	if( ( mm1 == 2 && dd1 == 29 && leapYear(yyyy1) || mm1 == 2 && dd1 == 28 && !leapYear(yyyy1) ) && 
		( mm2 == 2 && dd2 == 29 && leapYear(yyyy2) || mm2 == 2 && dd2 == 28 && !leapYear(yyyy2) ) )
		dd2 = 30;

	if( mm1 == 2 && dd1 == 29 && leapYear(yyyy1) || mm1 ==2 && dd1 == 28 && !leapYear(yyyy1) )
		dd1 = 30;

	if( dd2 == 31 && (dd1 == 30 || dd1 == 31) )
		dd2 = 30;

	if(dd1 == 31)
		dd1 = 30;

	long days = (yyyy2 - yyyy1) * 360 + (mm2 - mm1) * 30 + (dd2 - dd1);

	return days;
}


long cutility::dateOffset(const long dt)
{
	long yyyy = dt / 10000;
	long mm = (dt - yyyy * 10000) / 100;
	long dd = dt % 100;

	long mp, yp, T, offset;

	if(mm <= 2){
		mp = 0;
		yp = yyyy -1;
	}
	else{
		mp = static_cast<long>(0.4 * mm + 2.3);
		yp = yyyy;
	}
			
	T = (yp / 4) - (yp / 100) + (yp / 400);
	offset = 365 * yyyy + 31 * (mm -1) + dd + T - mp;

	return offset;
}


long cutility::numberOfCalendarDay( const long dt1, const long dt2) // dt is in format 'yyyymmdd'
{
	long off1 = dateOffset(dt1);
	long off2 = dateOffset(dt2);
	long days = gmax(off1, off2) - gmin(off1, off2);
		
	return days;
}

static double SIGN(double a, double b) 
{
	return b >= 0.0 ? fabs(a) : -fabs(a);
}

double cutility::zBrent(OBJ func, const double& l_, const double& h_, const double& epsi_)
{
	int i;
	double a = l_, b = h_, c = h_, d, e, min1, min2;
	double fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;
	const double EPS = 1.0e-16;
	const int ITMAX = 10000;

	if ((fa > epsi_ && fb > epsi_) || (fa < -epsi_ && fb < -epsi_)){
		std::cout << "Root must be bracketed in Brent method" << std::endl;
		return 0.25;
	}

	fc=fb;

	for (i = 0; i < ITMAX; ++i) {
		if ((fb > epsi_ && fc > epsi_) || (fb < -epsi_ && fc < -epsi_)) {
			c = a; //Rename a, b, c and adjust bounding interval d.
			fc = fa; 
			e = d = b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol1 = 2.0 * EPS * fabs(b) + 0.5 * epsi_; //Convergence check.
		xm = 0.5 * (c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) 
			return b;

		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb/fa;	//Attempt inverse quadratic interpolation.
			if (a == c){
				p = 2.0 * xm * s;
				q = 1.0 - s;
			} 
			else {
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (p > 0.0) 
				q = -q; //Check whether in bounds.
			p = fabs(p);
			min1 = 3.0 *xm * q - fabs(tol1 * q);
			min2 = fabs(e * q);
			if (2.0 * p < (min1 < min2 ? min1 : min2)) {
				e = d; //Accept interpolation.
				d = p / q;
			} 
			else{
				d = xm; //Interpolation failed, use bisection.
				e = d;
			}
		} 
		else{ //Bounds decreasing too slowly, use bisection.
			d = xm;
			e = d;
		}
		
		a = b; //Move last best guess to a.
		fa = fb;
		if (fabs(d) > tol1) //Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);
		fb = func(b);
	}
	
	std::cout << "Maximum number of iterations exceeded in Brent method" << std::endl;

	return 0.25; //Never get here.
}


double cutility::biSection(OBJ m_fcn, const double& l_, const double& h_, const double& epsi_)
// m_fcn must be an increasing function with m_fcn(l_) < 0 and m_fcn(h_) > 0
{
	double epsilon = epsi_;
	double low = l_, high = h_;
	
	double mid = (low+high)/2.0;

	double func = m_fcn(mid);

	while( fabs( func ) > epsilon){
		if( func > epsilon){
			high = mid;
		}
		else{
			low = mid;
		}

		mid = (low+high)/2.0;

		func = m_fcn(mid);

		if(high-low < epsilon)
			break;
	}

	return mid;
}


double cutility::normDensity(const double xp)
{
	const double root2pi_recip = 0.398942280401433; // == 1/sqrt(2*pi)
	double f_x = root2pi_recip * exp(-0.5 * xp * xp);

	return f_x;
}

// Computation N(z) = 1/sqrt(2*Pi) * int_{-inf}^{z}(exp(-t^2/2))
// error is around 1e-7
double cutility::normCDF(const double xp)
{
	double x = fabs(xp);

	const double a[] = {0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};	
	const double gamma = 0.2316419;
	const double root2pi_recip = 0.398942280401433; // == 1/sqrt(2*pi)

	double k = 1.0 /(1.0 + gamma * x);
	double f_x = root2pi_recip * exp(-0.5 * x * x);

	double z = 0.0;
	double kpow = k;
	for(int i=0; i<5; ++i){
		z += kpow * a[i];
		kpow *= k;
	}

	double v = 1.0 - z * f_x;

	if(xp < 0)
		v = 1.0 - v;

	return v;
}


double cutility::svCumNorm(const double xp)
{
	// Single varible normal CDF
	// From Applied Statistics algorithm 66, 
	// from the file normal.cpp of Matthew Willis
    double z = xp;
    double p,q;
    static double P [] = {
	220.2068679123761, 221.2135961699311, 112.0792914978709,
	33.91286607838300, 6.373962203531650, .7003830644436881,
	.03526249659989109 
    };
    static double Q [] = {  
	440.4137358247522, 793.8265125199484, 637.3336333788311,
	296.5642487796737, 86.78073220294608, 16.06417757920695, 
	1.755667163182642, .08838834764831844
    };
    static double CUTOFF = 7.071;
    static double root2pi = 2.506628274631001;
    double _pdf = 0.0;
    double zabs = z>0 ? z : -z;
    if (zabs >37.0) {
	if (z>0) {p = 1.; q = 0.;} else { q = 1.; p = 0.; }
	return p;
    }
    double expntl = exp(-0.5*zabs*zabs);
    _pdf = expntl/root2pi;
    if (zabs < CUTOFF) {
	p = expntl*((((((P[6]*zabs + P[5])*zabs + P[4])*zabs + P[3])*zabs +
		      P[2])*zabs + P[1])*zabs + P[0])/(((((((Q[7]*zabs + Q[6])*zabs +
                Q[5])*zabs + Q[4])*zabs + Q[3])*zabs + Q[2])*zabs + Q[1])*zabs +
                Q[0]);
    } else {
	p = _pdf/(zabs + 1.0/(zabs + 2.0/(zabs + 3.0/(zabs + 4.0/
						     (zabs + 0.65)))));
    }
    if (z<0) {
	q = 1-p;
    } else {
	q = p; p = 1-q;
    }
    return p;
}


double cutility::bivarCumNorm(const double a, const double b, const double r)
{
	// two dimension normal CDF

	int i;
	static double x[] = {
		0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992
	};

	static double w[] = {
		0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042
	};

	double h1, h2, LH, h12, h3, h5, h6, h7, h8;
	double r1, r2, r3, rr, AA, ab;

	double ret;

	h1 = a;
	h2 = b;
	h12 = (h1 * h1 + h2 * h2) / 2.0;

	LH = 0.0;

	if(fabs(r) >= 0.7){
		r2 = 1 - r * r;
		r3 = sqrt(r2);
		if( r < 0)
			h2 = -h2;
		h3 = h1 * h2;
		h7 = exp(-h3 / 2.0);
		if(fabs(r) < 1){
			h6 = fabs(h1 - h2);
			h5 = h6 * h6 / 2.0;
			h6 = h6 / r3;
			AA = 0.5 -h3 / 8;
			ab = 3.0 - 2.0 * AA * h5;
			LH = 0.13298076 * h6 * ab * (1 - svCumNorm(h6)) 
				- exp(-h5 / r2) * (ab + AA * r2) * 0.053051647;

			for(i=0; i<5; ++i){
				r1 = r3 * x[i];
				rr = r1 * r1;
				r2 = sqrt(1.0 - rr);
				if(h7 == 0)
					h8 = 0;
				else
					h8 = exp(-h3 / (1 + r2)) / (r2 * h7);
				LH = LH -w[i] * exp(-h5 / rr) * (h8 - 1 - AA * rr);
			}
		}
		
		ret = LH * r3 * h7 + svCumNorm(gmin(h1, h2));
		
		if(r < 0)
			ret = svCumNorm(h1) - ret;
	}
	else{
		h3 = h1 * h2;
		if(fabs(r) > 0){
			for(i=0; i<5; ++i){
				r1 = r * x[i];
				r2 = 1 - r1 * r1;
				LH = LH + w[i] * exp((r1 * h3 - h12) / r2) / sqrt(r2);
			}
		}

		ret = svCumNorm(h1) * svCumNorm(h2) + r * LH;
	}

	return ret;
}

double cutility::straddle(double tenor, double discount, double libor, 
						  double strike, double diffusion, double sigma)
// tenor is the distance from start to end in years
// discount is the discount factor at libor ends
// flibor is the forward libor rate used in the cap or floor 
// strike is the strike level
// diffusion is the distance from valuing date upto expiry in years (expiry may not be equal to the libor start)
// sigma is the volatility
{
	double d1 = (log(libor/strike) + 0.5 * sigma * sigma * diffusion)/(sigma * sqrt(diffusion));
	double d2 = d1 - sigma * sqrt(diffusion);

	double v = tenor * discount * (libor * (normCDF(d1) - normCDF(-d1)) + strike * (normCDF(-d2) - normCDF(d2)));

	return v;
}


void cutility::getRand_MB(double * randVec, int N)
// generate a vector of normal random variates, of size N, based on Marsaglia-Bray algorithm, 
// a modification of the Box-Muller algorithm.
{
	double sum;

	for (int i = 0; i < N; i += 2){
		double rand1, rand2;

		do{
			rand1 = 2.0 * (double)rand() / RAND_MAX - 1.0;
			rand2 = 2.0 * (double)rand() / RAND_MAX - 1.0;
			sum = rand1 * rand1 + rand2 * rand2;
		} while(sum >= 1);

		double c = sqrt( -2.0 * log(sum) / sum );
		randVec[i] = c * rand1;

		if(i < N-1) 
			randVec[i+1] = c * rand2;
	}
}


void cutility::getRand_Norm(double * randVec, int N)
// generate a vector of normal random variates, of size N, based on gasdev().
{
	double sum;

	for (int i = 0; i < N; ++i){
		double rand0 = (double)rand() / (RAND_MAX+1.0);
		randVec[i] = normalInverseCdf(rand0);
	}
}


double cutility::normalInverseCdf(double P)
{
    // From Applied Statistics PPND16, algorithm 241
    // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
    // Produces the normal deviate Z corresponding to a given lower
    // tail area of P; Z is accurate to about 1 part in 10**16.

	if(P >= 1.0 - 1e-16)
		P = 1.0 - 1e-16;

    static double zero = 0.0;
    static double one = 1.0;
    static double half = 0.5;
    static double split1=0.425;
    static double split2=5.0;
    static double const1=0.180625;
    static double const2=1.6;
    // Coefficients for P close to 0.5
    static double A[] = {	    
		3.3871328727963666080e0,	1.3314166789178437745e+2,
		1.9715909503065514427e+3,	1.3731693765509461125e+4,
		4.5921953931549871457e+4,	6.7265770927008700853e+4,
		3.3430575583588128105e+4,	2.5090809287301226727e+3 
	};
    static double B[] = {
		4.2313330701600911252e+1,	6.8718700749205790830e+2,
		5.3941960214247511077e+3,   2.1213794301586595867e+4,
		3.9307895800092710610e+4,	2.8729085735721942674e+4,
		5.2264952788528545610e+3
    };
    // Coefficients for P not close to 0, 0.5 or 1.
    static double C[] = { 
		1.42343711074968357734e0,	4.63033784615654529590e0,
		5.76949722146069140550e0,	3.64784832476320460504e0,
		1.27045825245236838258e0,	2.41780725177450611770e-1,
		2.27238449892691845833e-2,	7.74545014278341407640e-4 
    };    
    static double D[] = { 
		2.05319162663775882187e0,	1.67638483018380384940e0,
		6.89767334985100004550e-1,	1.48103976427480074590e-1,
		1.51986665636164571966e-2,	5.47593808499534494600e-4,
		1.05075007164441684324e-9
    };
    //	Coefficients for P near 0 or 1.
    static double E[] = { 
		6.65790464350110377720e0,    5.46378491116411436990e0,
		1.78482653991729133580e0,    2.96560571828504891230e-1,
		2.65321895265761230930e-2,   1.24266094738807843860e-3,
		2.71155556874348757815e-5,   2.01033439929228813265e-7
    };
    static double F[] = {
		5.99832206555887937690e-1,  1.36929880922735805310e-1,
		1.48753612908506148525e-2,  7.86869131145613259100e-4,
		1.84631831751005468180e-5,  1.42151175831644588870e-7,
		2.04426310338993978564e-15
    };
    
    double Q = P- half;
    double absQ = Q<0 ? -Q : Q;
    double R;
    double PPND16;

    if (absQ <= split1) {
		R = const1 - Q*Q;

		PPND16 = Q * (((((((A[7] * R + A[6]) * R + A[5]) * R + A[4]) * R + A[3])
     			* R + A[2]) * R + A[1]) * R + A[0]) /
		 (((((((B[6] * R + B[5]) * R + B[4]) * R + B[3]) * R + B[2])
			 * R + B[1]) * R + B[0]) * R + one);

		return PPND16;
    }

    if(Q < zero){       
		R = P;
    }
	else{
		R = one - P;
    }

    if (R <= zero) {
		return 0.0; // error
    }

    R = sqrt(-log(R));

    if(R <= split2){
		R = R - const2;	
		PPND16 = (((((((C[7] * R + C[6]) * R + C[5]) * R + C[4]) * R + C[3])
			    * R + C[2]) * R + C[1]) * R + C[0]) /
		 (((((((D[6] * R + D[5]) * R + D[4]) * R + D[3]) * R + D[2])
			   * R + D[1]) * R + D[0]) * R + one);
    } 
	else{
		R = R - split2;
		PPND16 = (((((((E[7] * R + E[6]) * R + E[5]) * R + E[4]) * R + E[3])
			    * R + E[2]) * R + E[1]) * R + E[0]) /
		 (((((((F[6] * R + F[5]) * R + F[4]) * R + F[3]) * R + F[2])
			  * R + F[1]) * R + F[0]) * R + one);
    }
	
    if(Q <= zero) 
		PPND16 = - PPND16;

    return PPND16;
}


