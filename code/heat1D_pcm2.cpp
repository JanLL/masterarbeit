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

//#include "ind_dyn_model_description.hpp"
//#include "ind_compile_time_info.hpp"
#include "ind_interface.hpp"
#include "ind_taylor_coeffs.hpp"

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


std::vector< Sonic::DMat> g_traj;
std::vector< Sonic::DMat> g_fwdSens;
std::vector< Sonic::DMat> g_adjSens;


const long maxLevels = 2; // Comment: for more levels, the boost preprocessor subset might be helpful

double L1;					// Length of Constantan
double L3;					// Length of PCM
long N1[maxLevels];			// Number of discretization points Constantan
long N3[maxLevels];			// Number of discretization points PCM
double dx_const[maxLevels];	// Mesh size Constantan
double dx_pcm[maxLevels];	// Mesh size PCM
double alpha;				// quotient for Constantan / PCM grid transition

template <typename T, long level>
svLong diffRHS(TArgs_ffcn<T> &args, TDependency *depends)
{

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



	const T* x = args.xd;

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


	// PCM
	T c_p_j;
	T dc_p_j;
	for (long j = N1[level]; j <= N1[level]+N3[level]-2; j++)
	{
		// c_p calculation using Fraser-Suzuki Peak
		//fraser_suzuki(x[j], c_p_j, dc_p_j, h, r, wr, sr, z, b);
	
		// c_p calculation using old atan formula
		c_p_formula(x[j], c_p_j, dc_p_j, args.p[4], args.p[5], args.p[6], args.p[7], args.p[8], args.p[9]);

		// linear part 
		args.rhs[j] = scale_pcm/c_p_j * (x[j-1] - 2.0 * x[j] + x[j+1]);

		// non-linear part (just c_p temperature dependent atm, lambda and rho constant)
		args.rhs[j] -= scale_pcm/(4.*c_p_j*c_p_j) * dc_p_j * (x[j+1] - x[j-1])*(x[j+1] - x[j-1]);
	}

	// RHS boundary, no flux
	args.rhs[N1[level]+N3[level]-1] = scale_pcm * (x[N1[level]+N3[level]-2] - x[N1[level]+N3[level]-1]); // Neumann boundary (right)


	return 0;
}





int main( int argc, char* argv[] )
{
	svLong errorCode ( 0 );


	// create a DAESOL-II instance
	IIntegrator * integrator = IIntegrator::create( IIntegrator::Impl_DAESOL_II );

	// get corresponding integrator options and evaluator from the integrator instance
	IOptions   *intOptions   = integrator-> getOptions();
	TIntegratorEvaluator *evaluator    = integrator-> getEvaluator();


	// choose the modules used by DAESOL-II by passing the corresponding options (ATLAS or UMFPACK for linear algebra)
	daesol2::CorrectorEqnSolverOptions* corrEqnSolvOpt = new daesol2::StdNewtonSolverATLASOptions();
	daesol2::BDFStarterOptions* bdfStartOpt = new daesol2::SelfStarterATLASOptions();
	//daesol2::CorrectorEqnSolverOptions* corrEqnSolvOpt = new daesol2::StdNewtonSolverUMFPACKOptions();
	//daesol2::BDFStarterOptions* bdfStartOpt = new daesol2::SelfStarterUMFPACKOptions();
	
	daesol2::ConsistencySolverOptions *consOpt = new daesol2::NewtonConsistencySolverOptions();


	intOptions->setOption(
			IOptions::Option_DAESOL_CorrectorEquationSolver,
			corrEqnSolvOpt
			);
	intOptions->setOption(
			IOptions::Option_DAESOL_BDFStarter,
			bdfStartOpt
			);
	intOptions->setOption(
		IOptions::Option_DAESOL_ConsistencySolver,
		consOpt
		);

	intOptions->setOption( IOptions::Option_DAESOL_INDMode, std::string("Iterative").c_str() );


	// create and configure dimensions
	TIntegratorDimensions dims;
	dims. dim[ Component_XD ] = 2;
	dims. dim[ Component_XA ] = 0;
	dims. dim[ Component_U  ] = 1;
	dims. dim[ Component_P  ] = 1;
	dims. dim[ Component_Q  ] = 1;
	dims. dim[ Component_H  ] = 1;
	dims. nTrajectories       = 1;



	return errorCode;
}
