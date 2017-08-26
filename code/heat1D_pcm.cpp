/**
 * \file heat1Drobin.cpp
 * \author Andreas Potschka
 * \date April 28, 2009
 *
 * \brief Finite Difference heat equation in 1D with Robin boundary control
 *
 */

#include <cmath>
#include <stdlib.h>
#include <complex>

#include "ind_dyn_model_description.hpp"
#include "ind_compile_time_info.hpp"
#include "sonic++.h"
#include <sstream>


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
	const T* x = args.xd;  // pointer to constant T (read only!)
	
	// parameters
	T a_const    = args.p[0];		// Temp.-conductivity (diffusion constant) of Constantan [mm^2/s] 
	T lambda_pcm = args.p[1];		// heat conductivity of pcm [mW/(mm*K)]
	T rho_pcm    = args.p[2];		// density of pcm [mg/mm^3]
	T heat_rate  = args.p[3];		// rate [K/min] with which oven temperature increases.

	T c_p_pcm = 2.; // [mJ/(mg*K)]   (testweise konstant erstmal)

	// some pre-calculations
	T heat_rate_s = heat_rate / 60; // [K/min] -> [K/s]
	T scale_Const = a_const / (dx_const[level]*dx_const[level]);
	T scale_pcm   = lambda_pcm / (rho_pcm*c_p_pcm);

	/******************* Building up RHS of ODE ********************/
	args.rhs[0] = heat_rate_s; // oven boundary with constant slope
	
	// Constantan linear part
	for (long j = 1; j <= N1[level]-1; j++)
	{
		args.rhs[j] = scale_Const * (x[j-1] - 2.0 * x[j] + x[j+1]); // Laplace 3-star
	}


	// intermediate area between Constantan and PCM
	int j = N1[level];
	args.rhs[N1[level]] = scale_Const * (2./(1.+alpha) * x[j-1] - 2./alpha * x[j] + 2./(alpha*(alpha+1.)) * x[j+1]);


	// PCM linear part
	for (long j = N1[level]+1; j < N1[level]+N3[level]-1; j++)
	{
		args.rhs[j] = scale_pcm * (x[j-1] - 2.0 * x[j] + x[j+1]); // Laplace 3-star
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
		test1 >> N1[level];  // for chosen level before, set grid size N1
		test1 >> N3[level];  // for chosen level before, set grid size N3
	}
	else  // standard choice of grid size
	{
		level = 0;
		N1[level] = 200;
		N3[level] = 50;
	}

	std::cout << "Using 1D heat equation with " << N1[level] <<" Discretization points for Constantan." << std::endl;


	dx_const[level] = L1 / static_cast<double>(N1[level]);
	dx_pcm[level] = L3 / static_cast<double>(N3[level]);
	alpha = L3/L1 * static_cast<double>(N1[level])/static_cast<double>(N3[level]);

	std::cout << "alpha: " << sprinf('%0.3f', alpha) << std::endl;

	m_dims. dim [ Component_T  ] = 1;
	m_dims. dim [ Component_XD ] = N1[level] + N3[level];
	m_dims. dim [ Component_XA ] = 0;
	m_dims. dim [ Component_P  ] = 4;
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


