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


#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/nurbs.hpp"
#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/nurbs.cpp"
#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/Interp1d_linear.hpp"
#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/Interp1d_linear.cpp"


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
size_t num_cntrl_pts;		// number of NURBS control points
const int nurbs_order = 4;	// NURBS curve will be in C^2 with nurbs_order=4

template <typename T, long level>
svLong heat_eq_rhs(TArgs_ffcn<T> &args, TDependency *depends)
{
	static std::vector<double> cntrl_pts_x(num_cntrl_pts, 0.);
	static std::vector<double> cntrl_pts_y(num_cntrl_pts, 0.);

	const T* x = args.xd;  // pointer to constant T (read only!)
	
	// parameters
	T a_const    = args.p[0];		// Temp.-conductivity (diffusion constant) of Constantan [mm^2/s] 
	T lambda_pcm = args.p[1];		// heat conductivity of pcm [mW/(mm*K)]
	T rho_pcm    = args.p[2];		// density of pcm [mg/mm^3]
	T heat_rate  = args.p[3];		// rate [K/min] with which oven temperature increases.


	T c_p_pcm = 2.; // [mJ/(mg*K)]   (testweise konstant erstmal)

	// Versuche den double value aus dem adouble raus zu bekommen...
	//T const * ptr = args.p;
	//(*ptr).value();
	// Hier gibts einen Konflikt weil value() weiterleitet zu getValue(), 
	// welches eine const methode ist und der compiler angst hat man wuerde
	// mit dem Pointer dann was aendern... obwohl es ein ptr auf konstantes adouble
	// ist. 
	// Ganz am Anfang der Errors sagt er, dass *ptr ein const double ist.

	//const T test = args.p[0];
	//test.getValue();
	// Hier kommt in der fehlermeldung, dass test vom typ const double sei 
	// (und daher keine methode .value() hat), aber es ist doch ein adouble....)

	//T const * ptr = args.p;
	//const double test = *ptr;
	// Hier sagt er dann allerdings, dass *ptr ein const adouble ist, im Gegensatz zum 1. Versuch

	/*std::ostringstream oss;
	oss << a_const;
	std::string s = oss.str();
	double a_const_double = std::stod(s, NULL);
	std::cout << a_const_double << std::endl;*/
	// Workaround der ganz und garnicht schoen ist, aber es funktioniert... ;D

	std::vector<double> cntrl_pts_x_input(num_cntrl_pts);
	std::vector<double> cntrl_pts_y_input(num_cntrl_pts);
	double c_p_param_input;
	std::ostringstream oss;
	std::string str;

	// Ist zwar echt nicht huebsch, aber laufzeittechnisch eig. kein Problem.
	// Pro Aufruf der RHS braucht die folgende for-Schleife 35us fuer 12 cntrl_pts.
	for (size_t i=0; i<num_cntrl_pts; ++i) {

		oss.str("");
		oss << args.p[4+i];
		str = oss.str();
		c_p_param_input = std::stod(str, NULL);
		cntrl_pts_x_input[i] = c_p_param_input;

		oss.str("");
		oss << args.p[4+num_cntrl_pts+i];
		str = oss.str();
		c_p_param_input = std::stod(str, NULL);
		cntrl_pts_y_input[i] = c_p_param_input;
	}


	// OLD: problem das man nicht direkt auf die parameter zugreifen kann...
	/*bool changed_cntrl_pts = false;
	for (size_t i=0; i<num_cntrl_pts; ++i) {
		if (cntrl_pts_x[i] != (args.p[4+i]).value() || 
			cntrl_pts_y[i] != (args.p[4+num_cntrl_pts+i]).value() ) {

			changed_cntrl_pts = true;
			break;
		}

	if (changed_cntrl_pts) {

		// Replace control points by new parameter values
		for (size_t i=0; i<num_cntrl_pts; ++i) {
			cntrl_pts_x[i] = (args.p[4+i]).value();
			cntrl_pts_y[i] = (args.p[4+num_cntrl_pts+i]).value();
		}

		// Build/Replace linear interpolation object for c_p(T) evaluation
		Nurbs nurbs(num_cntrl_pts, nurbs_order);
		nurbs.set_cntrl_pts_x(cntrl_pts_x);
		nurbs.set_cntrl_pts_y(cntrl_pts_y);

	}*/


	if (!(std::equal(cntrl_pts_x.begin(), cntrl_pts_x.end(), cntrl_pts_x_input.begin()) &&
	 	  std::equal(cntrl_pts_y.begin(), cntrl_pts_y.end(), cntrl_pts_y_input.begin()))) {
		// Bem.: Dieser if-Block wird bei Aenderung der cntrl_pts 2 mal aufgerufen -> unklar warum

		// Replace internal control point vectors
		cntrl_pts_x = cntrl_pts_x_input;
		cntrl_pts_y = cntrl_pts_y_input;

		// Build/Replace linear interpolation object for c_p(T) evaluation
		Nurbs nurbs(num_cntrl_pts, nurbs_order);
		nurbs.set_cntrl_pts_x(cntrl_pts_x);
		nurbs.set_cntrl_pts_y(cntrl_pts_y);

		double h = 0.001;
		std::vector<double> C_x(int(1/h)+1);
		std::vector<double> C_y(int(1/h)+1);

		std::vector<double> C_temp(2);
		int i=0;
		for (double u=0. + 1e-8; u <= 1.; u+=h) {
			C_temp = nurbs.eval_nurbs_curve(u);
			C_x[i] = C_temp[0];
			C_y[i] = C_temp[1];
		
			i++;
		}

		// 1D Interpolation part
		Interp1d_linear c_p_interpolator(C_x, C_y);

		// TEST interpolator
		//std::vector<double> c_p_i(2);
		//c_p_i = c_p_interpolator(130.);
		//std::cout << c_p_i[0] << "\t" << c_p_i[1] << std::endl;
	}


	// some pre-calculations
	T heat_rate_s = heat_rate / 60; // [K/min] -> [K/s]
	T scale_Const = a_const / (dx_const[level]*dx_const[level]);
	T scale_pcm   = lambda_pcm / (rho_pcm*c_p_pcm);

	std::vector<T> c_p_pcm(N3[level]);
	// TODO: c_p_pcm mit c_p_interpolator fuellen


	/******************* Building up RHS of ODE ********************/
	args.rhs[0] = heat_rate_s; // oven boundary with constant slope
	
	// Constantan linear part
	for (long j = 1; j <= N1[level]-1; j++)
	{
		args.rhs[j] = scale_Const * (x[j-1] - 2.0 * x[j] + x[j+1]); 
	}


	// intermediate area between Constantan and PCM
	int j = N1[level];
	args.rhs[N1[level]] = scale_Const * (2./(1.+alpha) * x[j-1] - 2./alpha * x[j] + 2./(alpha*(alpha+1.)) * x[j+1]);


	// PCM linear part
	for (long j = N1[level]+1; j < N1[level]+N3[level]-1; j++)
	{
		args.rhs[j] = scale_pcm * (x[j-1] - 2.0 * x[j] + x[j+1]); 
	}

	args.rhs[N1[level]+N3[level]-1] = scale_pcm * (x[N1[level]+N3[level]-2] - x[N1[level]+N3[level]-1]); // Neumann boundary (right)

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
		test1 >> num_cntrl_pts; // number of control points of c_p parametrization with NURBES.
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
	m_dims. dim [ Component_P  ] = 28;
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


