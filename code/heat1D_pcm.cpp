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
long N[maxLevels];	///< Number of discretization points
double h[maxLevels];	///< Mesh size

template <typename T, long level>
svLong heat_eq_rhs(TArgs_ffcn<T> &args, TDependency *depends)
{
	const T* x = args.xd;  // pointer to constant T (read only!)
	
	// parameters
	T L1         = args.p[0];	    // Length Constantan [mm]
	T L3         = args.p[1];		// Length PCM [mm]
	T a_const    = args.p[2];		// Temp.-conductivity (diffusion constant) of Constantan [mm^2/s] 
	T lambda_pcm = args.p[3];		// heat conductivity of pcm [mW/(mm*K)]
	T rho_pcm    = args.p[4];		// density of pcm [mg/mm^3]
	T heat_rate  = args.p[5];		// rate [K/min] with which oven temperature increases.



	T scale = 1. / (h[level]*h[level]);

	//args.rhs[0] = scale * (x[1] - x[0] + (beta*h[level])*(u - x[0])); // Robin boundary 
	args.rhs[0] = 0.1;
	for (long j = 1; j < N[level] - 1; j++)
	{
		args.rhs[j] = scale * (x[j-1] - 2.0 * x[j] + x[j+1]); // Laplace 3-star
	}
	args.rhs[N[level]-1] = scale * (x[N[level]-2] - x[N[level]-1]); // Neumann boundary (right)

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
		test1 >> N[level];  // set for chosen leven before grid size N
	}
	else
	{
		level = 0;
		N[level] = 10;
	}

	std::cout << "Using 1D heat equation with " << N[level] <<" discretization points." << std::endl;

	h[level] = 1.0 / static_cast<double>(N[level] - 1);
	m_dims. dim [ Component_T  ] = 1;
	m_dims. dim [ Component_XD ] = N[level];
	m_dims. dim [ Component_XA ] = 0;
	m_dims. dim [ Component_P  ] = 6;
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


