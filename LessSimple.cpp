#include <fstream>
#include <string>
#include "LessSimple.h"

bool isNumeric(string stringToCheck) // string or a number

{ 
      bool numeric = false;  

      if(stringToCheck.find_first_not_of("0123456789.-") == string::npos)   
            numeric = true;  

      return numeric;
}


//double linearInterp(const double& thisDay, const map<double, double>& crv) // linear interpolation
//
//{
//
//      map<double, double>::const_iterator pos, pos1, pos2;
//
// 
//
//      double m_amount = 0.0;
//
// 
//
//      if(crv.empty()){
//
//            return m_amount;
//
//      }
//
// 
//
//      if(thisDay <= crv.begin()->first){
//
//            m_amount = crv.begin()->second;
//
//      }
//
//      else if(thisDay >= crv.rbegin()->first){
//
//            m_amount = crv.rbegin()->second;
//
//      }
//
//      else{
//
//            pos=crv.find(thisDay);
//
//            if(pos != crv.end())
//
//                  m_amount = pos->second;
//
//            else{
//
//                  pos2 = crv.upper_bound(thisDay);
//
//                  pos1 = pos2;
//
//                  advance(pos1, -1);
//
// 
//
//                  if(pos1 != crv.end() && pos2 != crv.end()){
//
//                        double m_t1 = pos1->first;
//
//                        double m_t2 = pos2->first;
//
//                        double m_t = thisDay;
//
// 
//
//                        m_amount = pos1->second + (pos2->second - pos1->second) /(m_t2 - m_t1) * (m_t - m_t1);
//
//                  }
//
//            }
//
//      }
//
// 
//
//      return m_amount;
//
//}
//
// 
//
// 
//
double timeWeightedRate(const double& thisDay, const map<double, double>& crv) // time weighted linear interpolation

{

      map<double, double>::const_iterator pos, pos1, pos2;

 

      double m_rate = 0.0;

 

      if(crv.empty()){

            return m_rate;

      }

 

      if(thisDay <= crv.begin()->first){

            m_rate = crv.begin()->second;

      }

      else if(thisDay >= crv.rbegin()->first){

            m_rate = crv.rbegin()->second;

      }

      else{

            pos=crv.find(thisDay);

            if(pos != crv.end())

                  m_rate = pos->second;

            else{

                  pos2 = crv.upper_bound(thisDay);

                  pos1 = pos2;

                  advance(pos1, -1);

 

                  if(pos1 != crv.end() && pos2 != crv.end()){

                        double m_t1 = pos1->first;

                        double m_t2 = pos2->first;

                        double m_t = thisDay;

 

                        double fr = (m_t2 * pos2->second - m_t1 * pos1->second) /(m_t2 - m_t1);

 

                        m_rate = 1.0/m_t * (m_t1 * pos1->second + fr * (m_t - m_t1)); // time-weighted linear

                  }

            }

      }

 

      return m_rate;

}

double getForwardRate(const double& start, const double& end, const map<double, double>& crv) // forward rate
{
      double m_rate1 = timeWeightedRate(start, crv);
      if(fabs(start - end) < 1.0e-12)
            return m_rate1;

      double m_rate2 = timeWeightedRate(end, crv);
      double m_rate = (m_rate2 * end - m_rate1 * start) / (end - start);

      return m_rate;
}

double getDF(const double& thisDay, const map<double, double>& crv) // time weighted linear rates, and the discount factor
{
	  map<double, double>::const_iterator pos, pos1, pos2;

      double m_rate, m_df = 0.0;

      if(crv.empty())
           return m_df;

      if(thisDay <= crv.begin()->first){
            m_df = exp(- thisDay * crv.begin()->second);
      }

      else if(thisDay >= crv.rbegin()->first){
            m_df = exp(- thisDay * crv.rbegin()->second);
      }

      else{

            pos=crv.find(thisDay);

            if(pos != crv.end())

                  m_df = exp(- thisDay * pos->second);

            else{

                  pos2 = crv.upper_bound(thisDay);

                  pos1 = pos2;

                  advance(pos1, -1);

 

                  if(pos1 != crv.end() && pos2 != crv.end()){

                        double m_t1 = pos1->first;

                        double m_t2 = pos2->first;

                        double m_t = thisDay;

 

                        double fr = (m_t2 * pos2->second - m_t1 * pos1->second) /(m_t2 - m_t1);

 

                        m_rate = 1.0/m_t * (m_t1 * pos1->second + fr * (m_t - m_t1)); // time-weighted linear

                        m_df = exp(- m_rate * m_t);

                  }

            }

      }

 

      return m_df;

}


double timeWeightedVariance(const double& thisDay, const map<double, double>& vol) // time weighted linear interploation of variances
{
      map<double, double>::const_iterator pos, pos1, pos2;

      double m_var = 0.0;

      if(vol.empty()){
            return m_var;
      }

      if(thisDay <= vol.begin()->first){
            m_var = vol.begin()->second * vol.begin()->second;
      }

      else if(thisDay >= vol.rbegin()->first){
            m_var = vol.rbegin()->second * vol.rbegin()->second;
      }
      else{
            pos=vol.find(thisDay);
            if(pos != vol.end())
                  m_var = pos->second * pos->second;
            else{
                  pos2 = vol.upper_bound(thisDay);
                  pos1 = pos2;
                  advance(pos1, -1);

                  if(pos1 != vol.end() && pos2 != vol.end()){
                        double m_t1 = pos1->first;
                        double m_t2 = pos2->first;
                        double m_t = thisDay;
 
                        double fr = (m_t2 * pos2->second * pos2->second - m_t1 * pos1->second * pos1->second) /(m_t2 - m_t1);

                        m_var = 1.0/m_t * (m_t1 * pos1->second * pos1->second + fr * (m_t - m_t1)); // time-weighted linear
                  }
            }
	  }
	  
      return m_var;
}


double getForwardVar(const double& start, const double& end, const map<double, double>& vol) // forward variances
// we assume that the variances are interpolated by a time-weighted scheme
{
      double m_var1 = timeWeightedVariance(start, vol);

      if(fabs(start - end) < 1.0e-12)
            return m_var1;

      double m_var2 = timeWeightedVariance(end, vol);
      double m_var = (m_var2 * end - m_var1 * start) / (end - start);

      return cutility::gmax(m_var, 0.0);
}

// time weighted linear interpolation

double TimeWeightedLinearInterp(const map<double, double>& data, const double t)
{
	double temp;
	map<double, double>::const_iterator pos, pos1;

	if(t <= data.begin()->first)
		temp = data.begin()->second;
	else if(t >= data.rbegin()->first)
		temp = data.rbegin()->second;
	else
	{
		pos = data.find(t);
		if(pos != data.end())
			temp = pos->second;
		else
		{
			pos = data.upper_bound(t);
			pos1 = pos;
			--pos1;
			double slope = (pos->first * pos->second - pos1->first * pos1->second)/(pos->first - pos1->first);
			temp = (1/t)*(pos1->first * pos1->second + slope * (t - pos1->first));
		}
	}

	return temp;
}

// Dot Product

double DotProduct(const vector<double>& x, const vector<double>& y)
{
	double answer = 0;
	int xSize = x.size();
	int ySize = y.size();

	if (xSize != ySize)
		exit(EXIT_FAILURE);

	for (int i=0; i<xSize; ++i)
	{
		answer += x[i] * y[i];
	}

	return answer;
}

// Matrix Transpose

void MatrixTranspose(const vector<vector<double> >& A, vector<vector<double> >& B)
{
	int i, j;
	int M = A.size();
	int N = A[0].size();

	B.clear();
	B.resize(N);

	for (i=0; i<N; ++i)
	{
		B[i].resize(M);

		for (j=0; j<M; ++j)
			B[i][j] = A[j][i];

	}
}

// Matrix Multiplication

void MatrixMult(const vector<vector<double> >& A, const vector<vector<double> >& B, vector<vector<double> > & C)
{
	int i, j, k;
	int M = A.size();
	int N = B.size();
	int P = B[0].size();
	double sum;

	if (A[0].size() != B.size())
		exit(EXIT_FAILURE);
	
	C.clear();
	C.resize(M);
	for (i=0; i<M; ++i)
		C[i].resize(P);

	for (i=0; i<M; ++i)
	{
		for (j=0; j<P; ++j)
		{
			sum = 0;
			for (k=0; k<N; ++k)
				sum += A[i][k] * B[k][j];
			C[i][j] = sum;
		}
	}
}

// Cholesky Decomposition

void Cholesky(const vector<vector<double> >& corr, vector<vector<double> >& A)
{

	double sum;
	int N; 
	int i;
	int j;
	int k;
	
	N = corr.size();
	A.resize(N);
	for (i=0; i<N; i++)
		A[i].resize(N);

	// computes entries of matrix A
	for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++)
		{
			if (j<i)
			{
				sum = 0;
				for (k=0; k<j; ++k)
					sum += A[i][k] * A[j][k];
				A[i][j] = (corr[i][j] - sum)/(A[j][j]);
			}
			else if (j==i)
			{
				sum = 0;
				for (k=0; k<i; ++k)
					sum += A[i][k] * A[i][k];
				A[i][j] = sqrt(1 - sum);
			}
			else
				A[i][j] = 0;
		}
	}
}

// Rainbow Option

double Rainbow(const vector<double>& weight, const vector<map<double, double> >& rates, const  map<double, double>& discounts,
			const vector<map<double, double> >& vols, const vector<double>& fixings,
			const vector<double>& spots, const vector<vector<double> >& corr, double T, int M)
{

	int N = weight.size();
	vector<double> eps(N);
	vector<double> X(N);
	vector<double> termRates(N);
	vector<double> termVars(N);
	vector<double> scenarioPrice(N);	// for each Monte Carlo simulation

	vector<double> returns(N);
	vector<vector<double> > A;		// coefficient matrix
	double payoff = 0;
	double sum;
	double discount = getDF(T, discounts);
	int i;
	int j;
	int k;

	init_genrand(static_cast<unsigned long>(time(NULL)));

	for (i=0; i<N; ++i)
	{
		termRates[i] = TimeWeightedLinearInterp(rates[i], T);
		termVars[i] = timeWeightedVariance(T, vols[i]);
	}

	Cholesky(corr, A);

	for (k=0; k<M; ++k)			// M iterations
	{
		for (i=0; i<N; ++i)		// N underliers
		{
			X[i] = 0;
			eps[i] = gasdev();

			for (j=0; j<i+1; ++j)
				X[i] += eps[j] * A[i][j];
			
			scenarioPrice[i] = spots[i] * exp(X[i] * sqrt(T * termVars[i]) + T * termRates[i] - 0.5 * (T * termVars[i]));
			returns[i] = (scenarioPrice[i] - fixings[i])/fixings[i];
		}

		sort(returns.begin(), returns.end());
		sum = 0;
	
		for (i=0; i<N; ++i)
		{
			sum += weight[N-1-i] * returns[i];
		}
	
		payoff += max(sum, 0.0);
	}	
	
	return (payoff/M)*discount;
}

// Read Data (Rainbow)

bool GetData(vector<double>& weights, vector<map<double, double> >& rates, map<double, double>& discounts,
			 vector<map<double, double> >& vols, vector<double>& fixings,
			 vector<double>& spots, vector<vector<double> >& corr)
{
	weights.clear();
	spots.clear();
	fixings.clear();
	corr.clear();
	discounts.clear();
	rates.clear();
	vols.clear();

	string fileName = "C:\\C++ Projects\\Simple\\basket_info.txt";
	ifstream fin(fileName.c_str());
	if(!fin.is_open())
		return false;
	else
	{
		string fixing, spot, weight;
		while(fin >> fixing >> spot >> weight)
		{
			//cout << fixing << "\t" << spot << "\t" << weight << endl;
			if(!isNumeric(fixing))
				continue;

			fixings.push_back(atof(fixing.c_str()));
			spots.push_back(atof(spot.c_str()));
			weights.push_back(atof(weight.c_str()));			
		}	
	}

	fin.clear();
	fin.close();

	fileName = "C:\\C++ Projects\\Simple\\corrs.txt";	// corr
	fin.open(fileName.c_str());

	if(!fin.is_open())
		return false;
	else
	{
		double c1, c2, c3, c4, c5;
		while(fin >> c1 >> c2 >> c3 >> c4 >> c5)
		{
			vector<double> temp;
			temp.push_back(c1);
			temp.push_back(c2);
			temp.push_back(c3);
			temp.push_back(c4);
			temp.push_back(c5);
			corr.push_back(temp);
		}
	}

	fin.clear();
	fin.close();

	fileName = "C:\\C++ Projects\\Simple\\discount.txt";
	fin.open(fileName.c_str());

	if(!fin.is_open())
		return false;
	else
	{
		double m_date, m_rate;
		while(fin >> m_date >> m_rate)
			discounts.insert(make_pair(m_date, m_rate));
	}

	fin.clear();
	fin.close();

	/*vector<string> Assets;
	Assets.push_back();
	Assets.push_back();
	Assets.push_back();
	Assets.push_back();
	Assets.push_back();*/


	fileName = "C:\\C++ Projects\\Simple\\crvs.txt";
	fin.open(fileName.c_str());	

	if(!fin.is_open())
		return false;
	else
	{
		string time, rate;
		map<double, double> tMap;

		while (fin >> time >> rate)
		{
			if(!isNumeric(time))
			{
				if(!tMap.empty())
				{
					rates.push_back(tMap);
					tMap.clear();
				}
				continue;
			}

			tMap.insert(make_pair(atof(time.c_str()), atof(rate.c_str())));
		}
		rates.push_back(tMap);
	}

	fin.clear();
	fin.close();

	fileName = "C:\\C++ Projects\\Simple\\vols.txt";
	fin.open(fileName.c_str());	

	if(!fin.is_open())
		return false;
	else
	{
		map<double, double> tMap;
		string time, rate;

		while (fin >> time >> rate)
		{
			if(!isNumeric(time))
			{
				if(!tMap.empty())
				{
					vols.push_back(tMap);
					tMap.clear();
				}
				continue;
			}

			tMap.insert(make_pair(atof(time.c_str()), atof(rate.c_str())));
		}		
		vols.push_back(tMap);
	}

	return true;
}

// Asian Option

double Asian(int flag, double K, double T, const vector<map<double, double> >& rates, const vector<map<double, double> >& vols,
			 const map<double, double>& discounts, const vector<double>& weights, const vector<double>& fixings,
			const vector<double>& times, const vector<double>& spots, const vector<vector<double> >& corr, int N)
{
	
	int i, j, k, l;
	int m = times.size();
	int n = weights.size();

	double start, end, payoff, discount;

	vector<vector<double> > forRates(m);
	vector<vector<double> > forVars(m);
	vector<vector<double> > A;					// coefficient matrix

	vector<double> scenarioPrice(n);
	vector<double> eps(n);
	vector<double> X(n);

	for (j=0; j<m; ++j)
	{
		forRates[j].resize(n);
		forVars[j].resize(n);
	}

	init_genrand(static_cast<unsigned long>(time(NULL)));	// seed

	start = 0;
	for (j=0; j<m; ++j)		// computes forward rates and variances
	{
		end = times[j];
		for (i=0; i<n; ++i)
		{
			forRates[j][i] = getForwardRate(start, end, rates[j]);
			forVars[j][i] = getForwardVar(start, end, rates[j]);
		}
		start = end;
	}

	Cholesky(corr, A);
	payoff = 0;
	discount = getDF(T, discounts);

	for (k=0; k<N; ++k)					// N iterations
	{
		scenarioPrice = spots;
		start = 0;
		double outerSum = 0;

		for (j=0; j<m; ++j)				// m time steps
		{
			end = times[j];
			double dt = end - start;
			double innerSum = 0;

			for (i=0; i<n; ++i)			// n underliers
			{
				X[i] = 0;
				eps[i] = gasdev();

				for (l=0; l<i+1; ++l)
					X[i] += A[i][l] * eps[l];			// check this line
	
				scenarioPrice[i] *= exp(X[i] * sqrt(dt * forVars[j][i]) + dt * forRates[j][i] - 0.5 * (dt * forVars[j][i]));
				innerSum += weights[i] * (scenarioPrice[i]/fixings[i] - K);
			}
			outerSum += innerSum;
			start = end;
		}

		payoff += max(outerSum/m, 0.0);
	}

	return (payoff/N)*discount;
}