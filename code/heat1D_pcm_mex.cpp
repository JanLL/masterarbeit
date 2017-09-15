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


#include "mex.h"


//#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/nurbs.hpp"
//#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/nurbs.cpp"
//#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/Interp1d_linear.hpp"
//#include "/home/argo/masterarbeit/code/tools/nurbs_interpolation_c++/Interp1d_linear.cpp"



using namespace std;
using namespace SolvInd;


svLong INTEGRATOR_PRINT_LVL ( 0 );

std::vector<double> solGrid;
svULong nmp;

svULong nxd;
svULong np;

std::vector< Sonic::DMat> g_fwdSens;
std::vector< Sonic::DMat> g_adjSens;

Sonic::DMat g_traj2;




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






class TContFwdSensGetter : public IPlugin
{
public:

	TContFwdSensGetter()
	:
		IPlugin()
	{
		m_event_filter.reset();
		m_event_filter.set ( IPlugin::Event_OutputGrid );

		m_output_grid = solGrid;

		solGrid_idx = 0;

	}

	virtual svLong atEvent (
		const TEvent       event,
		const TOutputData *data
	)
	{
		//std::cout << "TContFwdSensGetter::atEvent called at time: " << data->m_time << std::endl;		

		const double* T_ptr = data->m_solution_xd[0];
		for (svULong i=0; i<nxd; ++i) {
			g_traj2(solGrid_idx, i) = *T_ptr;
			T_ptr++;
		}

		//std::cout << "TContFwdSensGetter: stored traj:\n" << g_traj[ g_traj.size() - 1 ];

		//g_fwdSens.push_back( Sonic::DMat ( Sonic::cDMat ( data->m_fwdSensitivities[0], data->m_fwdSensLeaDim, data->m_dims [ Component_XD ], m_nRays * m_fwdTCOrder) ) );

		//std::cout << "TContFwdSensGetter: stored fwdSens:\n" << g_fwdSens[ g_fwdSens.size() - 1 ];

		solGrid_idx++;
		return 0;
	}

	svULong m_fwdTCOrder;
	svULong m_nRays;

	svULong solGrid_idx;

};



double L1;				// Length of Constantan
double L3;				// Length of PCM
long N1;				// Number of discretization points Constantan
long N3;				// Number of discretization points PCM
double dx_Const;		// Mesh size Constantan
double dx_pcm;			// Mesh size PCM
double alpha;			// quotient for Constantan / PCM grid transition

double lambda_Const;
double rho_Const;
double c_p_Const;
double lambda_pcm;
double rho_pcm;

double m_pcm;			// mass of PCM sample [mg]

double heat_rate;
double T_0;
double T_end;



template <typename T>
svLong diffRHS(TArgs_ffcn<T> &args, TDependency *depends)
{


	const T* x = args.xd;

	// some pre-calculations
	T a_Const      = lambda_Const / (c_p_Const * rho_Const);
	T heat_rate_s  = heat_rate / 60; // [K/min] -> [K/s]
	T scale_Const  = a_Const / (dx_Const*dx_Const);
	T scale_pcm    = lambda_pcm / (rho_pcm * dx_pcm*dx_pcm);

	/******************* Building up RHS of ODE ********************/
	// Oven boundary with constant slope of Temp.
	args.rhs[0] = heat_rate_s;
	
	// Constantan (just linear part)
	for (long j = 1; j <= N1-2; j++)
	{
		args.rhs[j] = scale_Const * (x[j-1] - 2.0 * x[j] + x[j+1]); 
	}


	// Intermediate area between Constantan and PCM, belongs to Constantan
	long j = N1-1;
	args.rhs[j] = scale_Const * (2./(1.+alpha) * x[j-1] - 2./alpha * x[j] + 2./(alpha*(alpha+1.)) * x[j+1]);


	// PCM
	T c_p_j;
	T dc_p_j;
	for (long j = N1; j <= N1+N3-2; j++)
	{
		// c_p calculation using Fraser-Suzuki Peak
		fraser_suzuki(x[j], c_p_j, dc_p_j, args.p[0], args.p[1], args.p[2], args.p[3], args.p[4], args.p[5]);
	
		// c_p calculation using old atan formula
		//c_p_formula(x[j], c_p_j, dc_p_j, args.p[0], args.p[1], args.p[2], args.p[3], args.p[4], args.p[5]);

		// linear part 
		args.rhs[j] = scale_pcm/c_p_j * (x[j-1] - 2.0 * x[j] + x[j+1]);

		// non-linear part (just c_p temperature dependent atm, lambda and rho constant)
		args.rhs[j] -= scale_pcm/(4.*c_p_j*c_p_j) * dc_p_j * (x[j+1] - x[j-1])*(x[j+1] - x[j-1]);
	}

	// RHS boundary, no flux
	args.rhs[N1+N3-1] = scale_pcm * (x[N1+N3-2] - x[N1+N3-1]); // Neumann boundary (right)


	return 0;
}




svLong adjInjector( const svULong idx, svULong& adjLeaDim, const double*& adjTrajIn )
{
	//std::cout << "adjInjector: called with idx = " << idx << "\n";

	static Sonic::DMat injMat ( N1+N3, 2*solGrid.size() );
	injMat.zero();

	injMat(N1  ,2*idx  ) = 1;
	injMat(N1+1,2*idx+1) = 1;

	adjLeaDim = injMat.memDimRows();
	adjTrajIn = injMat.getFVector();

//	cout << "adjInjector: injection Matrix:\n" << injMat << endl;

	return 0;
} 





void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	svLong errorCode ( 0 );




	svULong gridsize = mxGetN(prhs[0]);
	const double* input_ptr = mxGetPr(prhs[0]);
	for (svULong i=0; i<gridsize; ++i) {
		solGrid.push_back(*input_ptr);
		input_ptr++;
	}
	nmp = solGrid.size();  // number measure points


	// Set Parameters
	L1 = 15.;
	L3 = 0.5;
	N1 = 200;
	N3 = 50;

	lambda_Const = 23.;		
	rho_Const    = 8.9;
	c_p_Const    = 0.41;

	lambda_pcm   = 0.96;	// heat conductivity of pcm [mW/(mm*K)]
	rho_pcm      = 0.85;	// density of pcm [mg/mm^3]

	m_pcm 	     = 10.47;   // mass of pcm [mg], here: sample 407

	heat_rate    = 10.;     // Oven temp heat rate [K/min]   
	T_0  	     = 10.;     // [degC]
	T_end  	     = 200.;    // [degC]

	

	g_traj2 = Sonic::DMat(nmp, N1+N3);


	// Some pre-calculations
	double heat_rate_s = heat_rate / 60.; // [K/min] -> [K/s]
	dx_Const = L1 / static_cast<double>(N1);
	dx_pcm = L3 / static_cast<double>(N3);
	alpha = L3/L1 * static_cast<double>(N1)/static_cast<double>(N3);

	const double t_0 = 0.;
	const double t_end = (T_end - T_0) / heat_rate_s;

	//for (double T=T_0; T < T_end; T += 0.01) {
	//	solGrid.push_back ( (T - T_0) / heat_rate_s );
	//}


	// create a DAESOL-II instance
	IIntegrator * integrator = IIntegrator::create( IIntegrator::Impl_DAESOL_II );

	// get corresponding integrator options and evaluator from the integrator instance
	IOptions             *intOptions = integrator-> getOptions();
	TIntegratorEvaluator *evaluator  = integrator-> getEvaluator();


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
	dims. dim[ Component_XD ] = N1+N3;
	dims. dim[ Component_XA ] = 0;
	dims. dim[ Component_U  ] = 0;
	dims. dim[ Component_P  ] = 6;
	dims. dim[ Component_Q  ] = 0;
	dims. dim[ Component_H  ] = 0;
	dims. nTrajectories       = 1;



	nxd = dims. dim [ Component_XD ];
	np  = dims. dim [ Component_P  ];


	// pass problem dimensions to evaluator and integrator
	evaluator->  setDimensions ( dims );
	integrator-> setDimensions ( dims );

	// create model function class and fill it with the rhs pointer and pass it to the evaluator
	TModelFunctions model;
	model. setFunction < Function_ffcn > ( diffRHS <double>  );
	model. setFunction < Function_ffcn > ( diffRHS <adouble>  );

	evaluator->setModelFunctions ( &model );
	integrator->setIntegratorPrintLevel( INTEGRATOR_PRINT_LVL );




	double* init_xd = new double [ nxd ];
	for (svULong i=0; i<nxd; ++i) {
		init_xd[i] = T_0;
	}

	double* params = new double [ np ];
	// old atan peak test parameters
	/*params[0] = 125.;
	params[1] = 10.;
	params[2] = 0.01;
	params[3] = 10.;
	params[4] = 0.003;
	params[5] = 2.;*/

	// fraser-suzuki Peak test parameters
	params[0] = 50.;
	params[1] = 15.;
	params[2] = 25.;
	params[3] = 0.2;
	params[4] = 125.;
	params[5] = 2.;


	integrator->setInitialValues(
			Component_XD,
			nxd, // leading dim of init_xd
			init_xd
			);
	integrator->setInitialValues(
			Component_P ,
			np,
			params
			);


	const double relTol ( 1e-06 );
	integrator->setStepsizeBounds    ( 0, 0   );
	integrator->setRelativeTolerance ( relTol );
	integrator->setInitialStepsize   ( 1e-06  );
	integrator->setMaxIntSteps       ( 1000   );
	integrator->setTimeHorizon       ( t_0, t_end   );

	// Adjoint Sensivity generation
	svULong nAdjDir = 2;
  	svULong nAdjDirTotal = solGrid.size() * nAdjDir;
	double* adjSensDir = new double [ nAdjDirTotal * nxd ]; // in direction N1+1 and N1+2
  	memset ( adjSensDir, 0, nAdjDirTotal * nxd * sizeof ( double ) );
 
  	// Forward Sensivity generation
	svULong nFwdDir = np;
  	double* fwdSensDir = new double [ (1 + nxd + np) * np ]; // in direction ( t,xd,(xa),q,p,h ), where the t part is in fact ignored
	memset ( fwdSensDir, 0, (1 + nxd + np) * np * sizeof ( double ) );

	for ( svULong i = 0; i < np; i++ ){
		fwdSensDir [ (1+nxd+np)*i + 1+nxd+i] = 1;
	}



	TContFwdSensGetter fwdSensGetter;
	fwdSensGetter.m_fwdTCOrder = 1;
	fwdSensGetter.m_nRays      = np;  // nicht sicher was das ist
	integrator->registerPlugin( &fwdSensGetter);
	

    integrator->activateFeature( IIntegrator::Feature_Store_Grid );
    integrator->activateFeature( IIntegrator::Feature_Store_Tape );  

	integrator->activateFeature( IIntegrator::Feature_Adjoint_Sensitivity_Injection );  
	integrator->setAdjointInjectionGrid( solGrid, &adjInjector );

    /*integrator->setForwardTaylorCoefficients (
			nFwdDir, // no of fwd Dirs
			1, // order
			1 + nxd + np, // leading Dim of fwd Dirs
			fwdSensDir
			);
	integrator->setAdjointTaylorCoefficients ( 0, 0, 0, 0, 0 );*/

	std::cout << "Integrate now...\n";
	errorCode = integrator->evaluate();
	std::cout << "Integrator return code : " << errorCode << "\n" << std::endl;
	if ( errorCode < 0 ){
		cout << "Error occured during evaluation, terminating now... \n" << errorCode<< std::endl;
		return;
	}
	

	integrator->setForwardTaylorCoefficients ( 0, 0, 0, 0 );
	integrator->setAdjointTaylorCoefficients (
			nAdjDirTotal, // no of adj Dirs
			0, // rays
			1, // order
			nxd, // leading Dim of adj Dirs
			adjSensDir
			);
  	std::cout << "Start backward sweep...\n";
	errorCode = integrator->backwardSensitivitySweep();
	std::cout << "Integrator return code : " << errorCode << "\n" << std::endl;
	if ( errorCode < 0 ){
		cout << "Error occured during beackward sweep, terminating now... \n" << errorCode<< std::endl;
		return;
	}


  	// Get Adjoint Sensitivities on solGrid and save in g_adjSens
	const double* sol;
	svULong solLD;

	for ( svULong ii = 0; ii < solGrid.size(); ii ++ ){
		g_adjSens.push_back( Sonic::DMat ( np, 2 ) );
		g_adjSens [ ii ].zero();
	}

	integrator->getAdjointSensitivities( 0, Component_P, solLD, sol );
	for ( svULong ii = 0; ii < solGrid.size(); ii ++ ){
		g_adjSens [ ii ] <<= Sonic::cDMat ( sol +  ii * 2 * solLD, solLD, np, 2 );
	}


	// Compute heat flux q between Constantan and PCM
	const double scale_q = (lambda_pcm * m_pcm) / (rho_pcm * dx_pcm*dx_pcm * N3); // Note: rho_pcm atm NOT temp.-dependend! Change Later!

	Sonic::DMat dq_dp (nmp, np);
	for (svULong i=0; i<nmp; ++i) {
		for (svULong j=0; j<np; ++j) {
			dq_dp(i,j) = scale_q * ( (g_adjSens[i])(j,0) - (g_adjSens[i])(j,1) );
		}
	}





	/******************** Output to Matlab *********************/
	// Temperature trajectory
	/*Sonic::DMat T_traj(nmp, N1+N3);
	for (svULong i=0; i<nmp; ++i) {
		T_traj.refRow(i) = g_traj[i].transpose();
	}
	const Sonic::cDMat& T_temp(T_traj);
	plhs[0] = mxCreateDoubleMatrix(T_temp.nRows(), T_temp.nCols(), mxREAL);
	Sonic::DMat T_output(mxGetPr(plhs[0]), T_temp.nRows(),
			T_temp.nRows(), T_temp.nCols(), Sonic::CreateReference);
	T_output <<= T_temp;*/

	const Sonic::cDMat& T_temp(g_traj2);
	plhs[0] = mxCreateDoubleMatrix(T_temp.nRows(), T_temp.nCols(), mxREAL);
	Sonic::DMat T_output(mxGetPr(plhs[0]), T_temp.nRows(),
			T_temp.nRows(), T_temp.nCols(), Sonic::CreateReference);
	T_output <<= T_temp;


	// Residuum Vector




	// Jacobian dq/dp
	const Sonic::cDMat& dq_dp_temp(dq_dp);
	plhs[1] = mxCreateDoubleMatrix(dq_dp_temp.nRows(), dq_dp_temp.nCols(), mxREAL);
	Sonic::DMat dq_dp_output(mxGetPr(plhs[1]), dq_dp_temp.nRows(),
			dq_dp_temp.nRows(), dq_dp_temp.nCols(), Sonic::CreateReference);
	dq_dp_output <<= dq_dp_temp;



	return;
}
