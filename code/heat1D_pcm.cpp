/**
 * \file heat1D_pcm.cpp
 * \author Jan Lammel
 * \date August 26, 2017
 *
 * \brief Finite Difference heat equation in 1D as DSC measurement simulation
 *
 */

#include <cmath>
#include <stdlib.h>
#include <complex>
#include <algorithm>
#include <string>
#include <time.h>

#include "ind_dyn_model_description.hpp"
#include "ind_compile_time_info.hpp"
#include "sonic++.h"
#include <sstream>
#include <fstream>


#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/nurbs.hpp"
#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/nurbs.cpp"
#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/Interp1d_linear.hpp"
#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/Interp1d_linear.cpp"


template<typename T>
void fraser_suzuki(T x, T& c_p, T& dc_p, T h,T r,T wr,T sr,T z, T b) {


	T c_p_case0;
	T c_p_case1;
	T dc_p_case0;
	T dc_p_case1;
	T condition;

	T log_sr = log(sr);

	// TODO: Mal mit fmin arbeiten statt komplett condition in condassign zu stecken...
	condition = z - wr*sr/(sr*sr-1) - x;

	T log_arg = (1 + (x-z)*(sr*sr-1)/(wr*sr));
	T exp_arg = -log(r)/(log_sr*log_sr) * log(log_arg)*log(log_arg);
	
	c_p_case0 = b;
	c_p_case1 = h*exp(exp_arg) + b;
	condassign(c_p, condition, c_p_case1, c_p_case0);

	dc_p_case0 = 0.;
	dc_p_case1 = -2.*log(r)/(log_sr*log_sr) * (sr*sr - 1)/(wr*sr) * log(log_arg)/log_arg * h*exp(exp_arg);
	condassign(dc_p, condition, dc_p_case1, dc_p_case0);

	return;
}


template<typename T>
void c_p_formula(T x, T& c_p, T& dc_p, T p0, T p1, T p2, T p3, T p4, T p5) {

	T atan_arg = -p3*(x-p0);
	T atan_value = atan(atan_arg);
	T exp_arg = -p2*(x-p0)*(x-p0);
	T exp_value = exp(exp_arg);

	c_p = (atan_value + M_PI/2.)*(p1*exp_value) + p4*x + p5;

	dc_p = -(p1*p3*exp_value)/(1+atan_value*atan_value) 
	       -2*p1*p2*(x-p0)*exp_value * (atan_value + M_PI/2.) + p4;

	return;
}




using namespace std;
using namespace SolvInd;

const long maxLevels = 2; // Comment: for more levels, the boost preprocessor subset might be helpful

double L1;					// Length of Constantan
double L3;					// Length of PCM
long N1[maxLevels];			// Number of discretization points Constantan
long N3[maxLevels];			// Number of discretization points PCM
double dx_const[maxLevels];	// Mesh size Constantan
double dx_pcm[maxLevels];	// Mesh size PCM
double alpha;				// quotient for Constantan / PCM grid transition

template <typename T, long level>
svLong heat_eq_rhs(TArgs_ffcn<T> &args, TDependency *depends)
{
	clock_t begin = clock();

	// parameters
	T a_const    = args.p[0];		// Temp.-conductivity (diffusion constant) of Constantan [mm^2/s] 
	T lambda_pcm = args.p[1];		// heat conductivity of pcm [mW/(mm*K)]
	T rho_pcm    = args.p[2];		// density of pcm [mg/mm^3]
	T heat_rate  = args.p[3];		// rate [K/min] with which oven temperature increases.

	T h          = args.p[4];
	T r          = args.p[5];
	T wr         = args.p[6];
	T sr         = args.p[7];
	T z          = args.p[8];
	T b          = args.p[9];

	/*double fs_params[6];
	std::ostringstream oss;
	oss << std::setprecision(20);
	std::string str;

	for (int i=0; i<=5; ++i) {
		oss.str("");
		oss << args.p[4+i];
		str = oss.str();
		fs_params[i] = std::stod(str, NULL);
	}*/

	// TEST c_p(T) and dc_p/dT(T) functions
	static bool one_time_call = true;
	if (one_time_call) {

		std::ofstream file_c_p;
  		file_c_p.open ("/home/argo/masterarbeit/code/c_p.txt");

  		T c_p;
  		T dc_p;
		for (T d=5.; d < 170; d += 0.01) {
			
			fraser_suzuki(d, c_p, dc_p, h, r, wr, sr, z, b);
			//c_p_formula(d, c_p, dc_p, args.p[4], args.p[5], args.p[6], args.p[7], args.p[8], args.p[9]);

			file_c_p << d << "\t" << c_p << "\t" << dc_p << std::endl; 

		}
		file_c_p.close();

		one_time_call = false;
	}


	const T* x = args.xd;  // pointer to constant T (read only!)
	



	// some pre-calculations
	T heat_rate_s  = heat_rate / 60; // [K/min] -> [K/s]
	T scale_Const  = a_const / (dx_const[level]*dx_const[level]);
	T scale_pcm    = lambda_pcm / (rho_pcm * dx_pcm[level]*dx_pcm[level]);

	/******************* Building up RHS of ODE ********************/
	// Oven boundary with constant slope of Temp.
	args.rhs[0] = heat_rate_s;
	
	// Constantan (just linear part)
	for (long j = 1; j <= N1[level]-2; j++)
	{
		args.rhs[j] = scale_Const * (x[j-1] - 2.0 * x[j] + x[j+1]); 
	}


	// Intermediate area between Constantan and PCM, belongs to Constantan
	long j = N1[level]-1;
	args.rhs[j] = scale_Const * (2./(1.+alpha) * x[j-1] - 2./alpha * x[j] + 2./(alpha*(alpha+1.)) * x[j+1]);


	//clock_t duration_fs = 0;
	// PCM
	T c_p_j;
	T dc_p_j;
	for (long j = N1[level]; j <= N1[level]+N3[level]-2; j++)
	{
		// c_p calculation using Fraser-Suzuki Peak

		/*oss.str("");
		oss << x[j];
		str = oss.str();
		double T_j = std::stod(str, NULL);
		double c_p_temp;
		double dc_p_temp;
		
		fraser_suzuki(T_j, c_p_temp, dc_p_temp, fs_params[0], fs_params[1], fs_params[2], fs_params[3], fs_params[4], fs_params[5]);
		c_p_j = c_p_temp;
		dc_p_j = dc_p_temp;

		std::cout << T_j << "\t" << c_p_j << "\t" << dc_p_j << std::endl;*/
		//clock_t begin_fs = clock();
		fraser_suzuki(x[j], c_p_j, dc_p_j, h, r, wr, sr, z, b);
		//clock_t end_fs = clock();
		//duration_fs += (end_fs - begin_fs);
		// c_p calculation using old atan formula
		//c_p_formula(x[j], c_p_j, dc_p_j, args.p[4], args.p[5], args.p[6], args.p[7], args.p[8], args.p[9]);

		// linear part 
		args.rhs[j] = scale_pcm/c_p_j * (x[j-1] - 2.0 * x[j] + x[j+1]);

		// non-linear part (just c_p temperature dependent atm, lambda and rho constant)
		args.rhs[j] -= scale_pcm/(4.*c_p_j*c_p_j) * dc_p_j * (x[j+1] - x[j-1])*(x[j+1] - x[j-1]);
	}

	// RHS boundary, no flux
	args.rhs[N1[level]+N3[level]-1] = scale_pcm * (x[N1[level]+N3[level]-2] - x[N1[level]+N3[level]-1]); // Neumann boundary (right)

	clock_t end = clock();
	//std::cout << std::setprecision(20) << double(end - begin) / CLOCKS_PER_SEC << std::endl;// "\t"
	          						   //<< double(duration_fs) / CLOCKS_PER_SEC << std::endl;

	return 0;
}

class myDynModel : public IDynamicModelDescription
{
public:
	myDynModel( const std::string& options = "");

private:
	static FFactory createMe;
	static IDynamicModelDescription::TRegisterTrigger myTrigger;
};

IDynamicModelDescription::TRegisterTrigger myDynModel::myTrigger( std::string("heat1D_pcm"), &myDynModel::createMe );
IDynamicModelDescription * myDynModel::createMe( const std::string& options )
{
	return new myDynModel( options );
}

myDynModel::myDynModel(  const std::string& options )
:
IDynamicModelDescription()
{
	long level = 0;

	// In diesem Block wird letztlich fuer ein Level im Mehrgitterverfahren die Anzahl der oertlichen Diskretisierungspunkte gewaehlt.
	string::size_type loc = options.find( "-", 0 );  // find location of char '-' in options-string of 'createDynamicModel' command
	if ( loc != string::npos ){  // check if char '-' exists at all
		std::istringstream test1 ( std::string ( options, 0, loc ) );  // load string stream with (here) numbers before the '-'
		test1 >> level;  // put first number in string stream into level
		if (level < 0 || level >= maxLevels)
		{
			std::cerr << "Error: Only " << maxLevels << " non-negative model levels available." << std::endl;
			return;
		}
		test1 >> L1;
		test1 >> L3;
		test1 >> N1[level];  	// for chosen level before, set grid size N1.
		test1 >> N3[level];     // for chosen level before, set grid size N3.
	}
	else  // standard choice of grid size
	{
		throw std::logic_error("Options string of importDynamicModelLib does not end with '-' !");
	}

	std::cout << "Using 1D heat equation with " << N1[level] <<" discretization points for Constantan." << std::endl;

	dx_const[level] = L1 / static_cast<double>(N1[level]);
	dx_pcm[level] = L3 / static_cast<double>(N3[level]);
	alpha = L3/L1 * static_cast<double>(N1[level])/static_cast<double>(N3[level]);

	std::cout << "alpha: " << alpha << std::endl;

	m_dims. dim [ Component_T  ] = 1;
	m_dims. dim [ Component_XD ] = N1[level] + N3[level];
	m_dims. dim [ Component_XA ] = 0;
	m_dims. dim [ Component_P  ] = 4 + 6;
	m_dims. dim [ Component_U  ] = 0;
	m_dims. dim [ Component_Q  ] = 0;
	//m_dims. dim [ Component_H  ] = 1;
	m_dims. nTrajectories        = 1;

	switch (level)
	{
		case 0:
			m_functions. setFunction<Function_ffcn>(heat_eq_rhs<double,0>);
			m_functions. setFunction<Function_ffcn>(heat_eq_rhs<adouble,0>);
			break;
		case 1:
			m_functions. setFunction<Function_ffcn>(heat_eq_rhs<double,1>);
			m_functions. setFunction<Function_ffcn>(heat_eq_rhs<adouble,1>);
			break;
	}
}


