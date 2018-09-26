// gutility.h

#ifndef GUTILITY_H_
#define GUTILITY_H_

#include <vector>

////////////// class for mean value ///////////////////////

class MeanValue{
private: 
  long num;
  double sum;
public:
  MeanValue() : num(0), sum(0.0){}

  void operator()(double elem){
    num++;
    sum += elem;
  }

  operator double(){
    return sum / static_cast<double>(num);
  }
};


typedef double (*OBJ)(double);

class cutility{
private: 
	static long dateOffset(const long dt);
	static int dayOfWeek(const long yyyymmdd);// 0 = sunday
	int mInt;
public:
	template<class T>
	static T gmin(T t1, T t2);

	template<class T>
	static T gmax(T t1, T t2);

	static bool leapYear(const long yyyy);

	static bool isWeekDay(const long yyyymmdd);

	static bool isWeekend(const long yyyymmdd);

	static long days30360( const long dt1, const long dt2);
	
	static long numberOfCalendarDay( const long dt1, const long dt2); // dt is in format 'yyyymmdd'

	static double biSection(OBJ m_fcn, const double& l_, const double& h_, const double& epsi_);

	static double zBrent(OBJ func, const double& l_, const double& h_, const double& epsi_);

	static double normDensity(const double xp);

	static double normCDF(const double xp);

	static double svCumNorm(const double xp); // single variable cumulative normal distribution

	static double bivarCumNorm(const double x, const double y, const double r); //bi-variate cum

	static double straddle(double tenor, double discount, double libor, 
							double strike, double diffusion, double sigma);
	static void getRand_MB(double * randVec, int N);

	static double normalInverseCdf(double prob);

	static void getRand_Norm(double * randVec, int N);

};


template<class T>
T cutility::gmin(T t1, T t2)
{
	return (t1>t2)? t2 : t1;
}

template<class T>
T cutility::gmax(T t1, T t2)
{
	return (t1>t2)? t1 : t2;
}


#endif
