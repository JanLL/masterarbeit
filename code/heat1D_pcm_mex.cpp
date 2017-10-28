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
#include "/home/argo/masterarbeit/code/c_p_parametrizations.cpp"



using namespace std;
using namespace SolvInd;

svLong INTEGRATOR_PRINT_LVL ( 0 );


enum c_p_param_types {

	uninitialized	  = 0,
	old_atan_formula  = 1,
	fraser_suzuki     = 2,
	gauss_linear_comb = 3,
	NURBS 			  = 4

};
c_p_param_types c_p_param_type;




// Global SolvIND objects
IIntegrator*          integrator;
TIntegratorEvaluator* evaluator;
TModelFunctions       model;


svULong nAdjDir;
svULong nAdjDirTotal;
double* adjSensDir;

svULong nFwdDir;
double* fwdSensDir;

std::vector<double> solGrid;
Sonic::DMat* q_meas;

svULong nmp;
svULong nxd;
svULong np;

Sonic::DMat* g_traj;
std::vector< Sonic::DMat> g_fwdSens;
std::vector< Sonic::DMat> g_adjSens;

bool initialized = false;


// Global Simulation parameters
double L1;					// Length of Constantan
double L3;					// Length of PCM
svULong N1;					// Number of discretization points Constantan
svULong N3;					// Number of discretization points PCM

/****** old **********/
double dx_Const;			// Mesh size Constantan
double dx_pcm;				// Mesh size PCM
double alpha;				// quotient for Constantan / PCM grid transition

/****** new **********/
// Spatial discretization gridsize vector of size [N1+N3-1]
std::vector<double> spatial_gridsize;
// vector of consecutive spatial discretization gridsize ratios for 2-point formula of 2nd derivative.
std::vector<double> alpha_vec;				

double lambda_Const;		// heat conductivity of Constantan [mW/(mm*K)]
double rho_Const;			// density of Constantan [mg/mm^3]
double c_p_Const;			// specific heat capacity of Constantan [mJ/(mg*K)]
double lambda_pcm;			// heat conductivity of PCM [mW/(mm*K)]
double rho_pcm;				// density of PCM [mg/mm^3]

double m_pcm;				// mass of PCM sample [mg]

double heat_rate;			// Oven heat rate [K/min]
double T_0;					// Start temperature of Constantan & PCM
double T_end;				// Integration ends when oven reaches T_end

// Global optimization parameter
double* c_p_params;



class TContFwdSensGetter : public IPlugin
{
public:
	TContFwdSensGetter()
	:
		IPlugin()
	{
		m_event_filter.reset();
		m_event_filter.set ( IPlugin::Event_OutputGrid );

		solGrid_idx = 0;

	}

	void setGrid(std::vector<double>& grid) {

		m_output_grid = grid;

	}

	void resetSolGrid_idx() {
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
			(*g_traj)(solGrid_idx, i) = *T_ptr;
			T_ptr++;
		}

		//std::cout << "TContFwdSensGetter: stored traj:\n" << g_traj[ g_traj.size() - 1 ];

		g_fwdSens.push_back( Sonic::DMat ( Sonic::cDMat ( data->m_fwdSensitivities[0], data->m_fwdSensLeaDim, data->m_dims [ Component_XD ], m_nRays * m_fwdTCOrder) ) );

		//std::cout << "TContFwdSensGetter: stored fwdSens:\n" << g_fwdSens[ g_fwdSens.size() - 1 ];

		solGrid_idx++;
		return 0;
	}

	svULong m_fwdTCOrder;
	svULong m_nRays;

	svULong solGrid_idx;

};

TContFwdSensGetter fwdSensGetter;  // aside other global variables because class declaration of TContFwdSensGetter needed first




svLong adjInjector( const svULong idx, svULong& adjLeaDim, const double*& adjTrajIn )
{
	//std::cout << "adjInjector: called with idx = " << idx << "\n";

	static Sonic::DMat injMat ( N1+N3, 2*solGrid.size() );
	static svULong N1_static = N1;
	static svULong N3_static = N3;
	static svULong solGridSize_static = solGrid.size();

	if (!(N1_static == N1 && N3_static == N3 && solGridSize_static == solGrid.size())) {
		injMat = Sonic::DMat( N1+N3, 2*solGrid.size() );
		N1_static = N1;
		N3_static = N3;
		solGridSize_static = solGrid.size();
	}
	
	injMat.zero();

	injMat(N1  ,2*idx  ) = 1;
	injMat(N1+1,2*idx+1) = 1;

	adjLeaDim = injMat.memDimRows();
	adjTrajIn = injMat.getFVector();

//	cout << "adjInjector: injection Matrix:\n" << injMat << endl;

	return 0;
} 


template <typename T>
void rho_pcm_formula(T x, T& rho, T& drho) {
// Formula and parameter from PCM_rho.m from Robert

	T p0 = 167.2182;
	T p1 = 129.9861;
	T p2 = 760.8218;
	T p3 = 0.078916;

	T exp_val = exp(p3*(x-p1));

	rho = p0 / (1 + exp_val) + p2;
	rho *= 1e-3;  // rescale to [mg/mm^3]

	drho = -p0*p3*exp_val / ((1+exp_val)*(1+exp_val));
	drho *= 1e-3; // rescale to [mg/mm^3]

	return;

}




template <typename T>
svLong diffRHS(TArgs_ffcn<T> &args, TDependency *depends)
{

	const T* x = args.xd;

	// some pre-processing
	T a_Const      = lambda_Const / (c_p_Const * rho_Const);
	T heat_rate_s  = heat_rate / 60; // [K/min] -> [K/s]
	T scale_Const  = a_Const / (dx_Const*dx_Const);
	T scale_pcm    = lambda_pcm / (dx_pcm*dx_pcm);


	/******************* Building up RHS of ODE ********************/
	// Oven boundary with constant slope of Temp.
	args.rhs[0] = heat_rate_s;
	
	// Constantan (just linear part)
	for (svULong j = 1; j <= N1-2; j++)
	{
		args.rhs[j] = scale_Const * (x[j-1] - 2.0 * x[j] + x[j+1]); 
	}


	// Intermediate area between Constantan and PCM, belongs to Constantan
	svULong j = N1-1;
	args.rhs[j] = scale_Const * (2./(1.+alpha) * x[j-1] - 2./alpha * x[j] + 2./(alpha*(alpha+1.)) * x[j+1]);


	// Set function pointer to chosen c_p parametrization function
	void (*c_p_fct)(T, T&, T&, const T*);
	switch (c_p_param_type)
		{
			case old_atan_formula:
				c_p_fct = &c_p_formula;
				break;
			case fraser_suzuki: 
				c_p_fct = &fraser_suzuki_formula;
				break;
			case gauss_linear_comb:
				c_p_fct = &gauss_linear_comb_formula;
				break;

			case NURBS:
				// Note: ADOL-C Index problem...
				mexErrMsgTxt( "NURBS c_p parametrization not implemented yet!" );
				break;
			default:
				mexErrMsgTxt( "Wrong c_p parametrization type!" );
		}


	// PCM
	T rho_j;
	T drho_j;
	T c_p_j;
	T h_j;
	for (svULong j = N1; j <= N1+N3-2; j++)
	{
		c_p_fct(x[j], c_p_j, h_j, args.p);
		
		/*if (j == N1) {
			std::cout << x[j] << "\t" << c_p_j << "\t" << h_j << "\n";
		}*/

		// compute density of pcm
		rho_pcm_formula(x[j], rho_j, drho_j);
		rho_j = 0.85; //  [mg/mm^3]

		// linear part 
		args.rhs[j] = scale_pcm/(c_p_j*rho_j) * (x[j-1] - 2.0 * x[j] + x[j+1]);
		//args.rhs[j] = scale_pcm * (1./(c_p_j*rho_j) + 1./(h_j*drho_j)) * (x[j-1] - 2.0 * x[j] + x[j+1]);

		// non-linear part (just c_p temperature dependent atm, lambda and rho constant) WRONG!
		//args.rhs[j] -= scale_pcm/(4.*c_p_j*c_p_j * rho_j) * dc_p_j * (x[j+1] - x[j-1])*(x[j+1] - x[j-1]);
		//args.rhs[j] -= scale_pcm/(4.*rho_j*rho_j * c_p_j) * drho_j * (x[j+1] - x[j-1])*(x[j+1] - x[j-1]);
		
	}

	// RHS boundary, no flux
	c_p_fct(x[N1+N3-2], c_p_j, h_j, args.p);
	//rho_pcm_formula(x[N1+N3-2], rho_j, drho_j);
	rho_j = 0.85; //  [mg/mm^3]

	args.rhs[N1+N3-1] = scale_pcm/(c_p_j*rho_j) * (x[N1+N3-2] - x[N1+N3-1]); // Neumann boundary (right)
	//args.rhs[N1+N3-1] = scale_pcm * (1./(c_p_j*rho_j) + 1./(h_j*drho_j)) * (x[N1+N3-2] - x[N1+N3-1]); // Neumann boundary (right)


	return 0;
}






void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	/****************************************************
	Function called by Matlab to run measurement simulation to get
	heat flux residuum vector and Jacobian w.r.t. c_p parameters.

	INPUT:
		\param prhs[0]: string, function call mode, 
			possible values: ["init", "evaluate", "optimization"]

		\param prhs[0]: measurement grid where ODE is evaluated.
		\param prhs[1]: simulation 

	OUTPUT:
		\param plhs[0]: 

	*****************************************************/ 

	char command[100];
	svLong errorCode ( 0 );
	if (nrhs < 1)
	{
		mexErrMsgTxt("heat1D_pcm cannot be called without command string as first argument!");
		return;
	}


	mxGetString(prhs[0], command, 99);

	// For Testing purposes...
	if (strncmp(command, "test", 99) == 0) {

	}




	if (strncmp(command, "init", 99) == 0) {

		if (initialized == true) {
			mexErrMsgTxt("Initialization was already done! You need to \"reset\" before reapply \"init\"!");
		}
		else {
			std::cout << "Initializing..." << std::endl;
		}

		// Simulation Parameters
		int n_sim_params = mxGetN(prhs[1]);
		if (n_sim_params != 12) {
			mexErrMsgTxt("Simulation parameter number inconstistent in C++ and Matlab!");
		}
		const double* input_ptr = mxGetPr(prhs[1]);
		L1 = *input_ptr; input_ptr++;
		L3 = *input_ptr; input_ptr++;
		N1 = *input_ptr; input_ptr++;
		N3 = *input_ptr; input_ptr++;

		lambda_Const = *input_ptr; input_ptr++;
		rho_Const    = *input_ptr; input_ptr++;
		c_p_Const    = *input_ptr; input_ptr++;

		lambda_pcm   = *input_ptr; input_ptr++;	  	// heat conductivity of pcm [mW/(mm*K)]
		m_pcm 	     = *input_ptr; input_ptr++;     // mass of pcm [mg], here: sample 407

		heat_rate    = *input_ptr; input_ptr++;     // Oven temp heat rate [K/min]   
		T_0  	     = *input_ptr; input_ptr++;     // [degC]
		T_end  	     = *input_ptr; input_ptr++;     // [degC]


		if (mxGetM(prhs[2]) != N1+N3-1) {
			mexErrMsgTxt("Spatial discretization grid size inconstistent with N1 and N3!");
		}
		input_ptr = mxGetPr(prhs[2]);
		for (svULong i=0; i<N1+N3-1; ++i) {
			spatial_gridsize.push_back(*input_ptr);
			input_ptr++;
		}


		// Measurement time grid and corresponding heat flux values
		nmp = mxGetM(prhs[3]);
		input_ptr = mxGetPr(prhs[3]);
		// Time grid of measurements (= for value and sensitivity extraction of ODE solution)
		for (svULong i=0; i<nmp; ++i) {
			solGrid.push_back(*input_ptr);
			input_ptr++;
		}
		// Heat flux measurement values
		q_meas = new Sonic::DMat(nmp,1);
		for (svULong i=0; i<nmp; ++i) {
			(*q_meas)(i,0) = *input_ptr;
			input_ptr++;
		}

		mxGetString(prhs[4], command, 99);
		if (strncmp(command, "old_atan_formula", 99) == 0) {
			c_p_param_type = old_atan_formula;
			np = 6;
		} else if (strncmp(command, "fraser_suzuki", 99) == 0) {
			c_p_param_type = fraser_suzuki;
			np = 6;
		} else if (strncmp(command, "gauss_linear_comb", 99) == 0) {
			c_p_param_type = gauss_linear_comb;
			np = 3*10 + 2 + 1;  // 10 Gaussians, linear, constant, start_enthalpy

		} else if (strncmp(command, "NURBS", 99) == 0) {
			mexErrMsgTxt( "NURBS problem with index unsolved..." ); // TODO!

			/*num_cntrl_pts = mxGetN(prhs[4]);
			cntrl_pts_x.resize(num_cntrl_pts);
			const double* input_ptr = mxGetPr(prhs[4]);
			for (svULong i=0; i<num_cntrl_pts; ++i) {
				cntrl_pts_x[i] = *input_ptr;
				input_ptr++;
			}
			np = num_cntrl_pts;

			
			c_p_param_type = NURBS;*/
		}
		c_p_params = new double [ np ];

		




		// Some pre-calculations / auxiliary variables
		g_traj = new Sonic::DMat(nmp, N1+N3);

		double heat_rate_s = heat_rate / 60.; // [K/min] -> [K/s]
		dx_Const = L1 / static_cast<double>(N1);
		dx_pcm = L3 / static_cast<double>(N3);
		alpha = L3/L1 * static_cast<double>(N1)/static_cast<double>(N3);

		const double t_0 = 0.;
		const double t_end = (T_end - T_0) / heat_rate_s;

		// create a DAESOL-II instance
		integrator = IIntegrator::create( IIntegrator::Impl_DAESOL_II );

		// get corresponding integrator options and evaluator from the integrator instance
		IOptions             *intOptions = integrator-> getOptions();
		evaluator 				         = integrator-> getEvaluator();


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
		dims. dim[ Component_P  ] = np;
		dims. dim[ Component_Q  ] = 0;
		dims. dim[ Component_H  ] = 0;
		dims. nTrajectories       = 1;


		nxd = dims. dim [ Component_XD ];


		// pass problem dimensions to evaluator and integrator
		evaluator->  setDimensions ( dims );
		integrator-> setDimensions ( dims );

		// create model function class and fill it with the rhs pointer and pass it to the evaluator
		model. setFunction < Function_ffcn > ( diffRHS <double>  );
		model. setFunction < Function_ffcn > ( diffRHS <adouble>  );

		evaluator->setModelFunctions ( &model );
		integrator->setIntegratorPrintLevel( INTEGRATOR_PRINT_LVL );




		double* init_xd = new double [ nxd ];
		for (svULong i=0; i<nxd; ++i) {
			init_xd[i] = T_0;
		}


		integrator->setInitialValues(
				Component_XD,
				nxd, // leading dim of init_xd
				init_xd
				);


		const double relTol ( 1e-7 );
		integrator->setStepsizeBounds    ( 0, 0   );
		integrator->setRelativeTolerance ( relTol );
		integrator->setInitialStepsize   ( 1e-06  );
		integrator->setMaxIntSteps       ( 1000   );
		integrator->setTimeHorizon       ( t_0, t_end   );

		// Adjoint Sensivity generation
		nAdjDir = 2;
	  	nAdjDirTotal = solGrid.size() * nAdjDir;
		adjSensDir = new double [ nAdjDirTotal * nxd ]; // in direction N1+1 and N1+2
	  	memset ( adjSensDir, 0, nAdjDirTotal * nxd * sizeof ( double ) );
	 
	  	// Forward Sensivity generation
		nFwdDir = np;
	  	fwdSensDir = new double [ (1 + nxd + np) * np ]; // in direction ( t,xd,(xa),q,p,h ), where the t part is in fact ignored
		memset ( fwdSensDir, 0, (1 + nxd + np) * np * sizeof ( double ) );

		for ( svULong i = 0; i < np; i++ ){
			fwdSensDir [ (1+nxd+np)*i + 1+nxd+i] = 1;
		}


		fwdSensGetter.setGrid(solGrid);
		fwdSensGetter.m_fwdTCOrder = 1;
		fwdSensGetter.m_nRays      = np;
		integrator->registerPlugin( &fwdSensGetter);
		

	    integrator->activateFeature( IIntegrator::Feature_Store_Grid );
	    integrator->activateFeature( IIntegrator::Feature_Store_Tape );  

		integrator->activateFeature( IIntegrator::Feature_Adjoint_Sensitivity_Injection );  
		integrator->setAdjointInjectionGrid( solGrid, &adjInjector );

		initialized = true;


	} else if (strncmp(command, "optimization", 99) == 0) {

		if (initialized == false) {
			mexErrMsgTxt("Not yet initialized. Use first command \"init\" before perform optimization!");
		}

		if (np != mxGetN(prhs[1])) {
			mexErrMsgTxt(" Wrong number of c_p input parameters! ");
		}
		const double* input_ptr = mxGetPr(prhs[1]);
		// Time grid of measurements (= for value and sensitivity extraction of ODE solution)
		for (svULong i=0; i<np; ++i) {
			c_p_params[i] = *input_ptr;
			input_ptr++;
		}


		evaluator->resetStatistics();
		fwdSensGetter.resetSolGrid_idx();
		g_adjSens.clear();
		g_fwdSens.clear();

		integrator->setAdjointInjectionGrid(solGrid, &adjInjector); 
		// Q: nicht ganz klar warum man das nochmal braucht? Sonst sind alle adjSens gleich 0...



		integrator->setInitialValues(Component_P, np, c_p_params);


	    integrator->setForwardTaylorCoefficients (
				nFwdDir,      // no of fwd Dirs
				1,            // order
				1 + nxd + np, // leading Dim of fwd Dirs
				fwdSensDir    // (1+nxd+np)x(np) matrix with np fwdDirections
				);
		//integrator->setForwardTaylorCoefficients ( 0, 0, 0, 0 );
		integrator->setAdjointTaylorCoefficients ( 0, 0, 0, 0, 0 );
		std::cout << "Integrate now...\n";
		clock_t t_begin = clock();
		errorCode = integrator->evaluate();
		clock_t t_end = clock();
		clock_t t_duration = t_end - t_begin;
		std::cout << "Integrator return code : " << errorCode << std::endl;
		std::cout << "Integration took " << double(t_duration) / CLOCKS_PER_SEC << " seconds." << std::endl;
		if ( errorCode < 0 ){
			cout << "Error occured during evaluation, terminating now... \n" << errorCode<< std::endl;
			return;
		}




		/*integrator->setForwardTaylorCoefficients ( 0, 0, 0, 0 );
		integrator->setAdjointTaylorCoefficients (
				nAdjDirTotal, // no of adj Dirs
				0, // rays
				1, // order
				nxd, // leading Dim of adj Dirs
				adjSensDir
				);
	  	std::cout << "Start backward sweep...\n";
		t_begin = clock();
		errorCode = integrator->backwardSensitivitySweep();
		t_end = clock();
		t_duration = t_end - t_begin;
		std::cout << "Integrator return code : " << errorCode << std::endl;
		std::cout << "Backward Sweep took " << double(t_duration) / CLOCKS_PER_SEC << " seconds." << std::endl;
		if ( errorCode < 0 ){
			cout << "Error occured during beackward sweep, terminating now... \n" << errorCode<< std::endl;
			return;
		}

	  	// Get Adjoint Sensitivities on solGrid and save in g_adjSens
		const double* sol;
		svULong solLD;

		for ( svULong ii = 0; ii < solGrid.size(); ii ++ ){
			g_adjSens.push_back( Sonic::DMat ( np, nAdjDir ) );
			g_adjSens [ ii ].zero();
		}

		integrator->getAdjointSensitivities( 0, Component_P, solLD, sol );
		for ( svULong ii = 0; ii < solGrid.size(); ii ++ ){
			g_adjSens [ ii ] <<= Sonic::cDMat ( sol +  ii * nAdjDir * solLD, solLD, np, nAdjDir );
		}

		// Computation of dq/dp = dq/dT * dT/dp with adjoint Sensitivities
		Sonic::DMat dq_dp_adj (nmp, np);
		for (svULong i=0; i<nmp; ++i) {
			
			double rho_pcm;
			double drho_pcm;

			double T_N1_i = (*g_traj)(i, N1);
			rho_pcm_formula(T_N1_i, rho_pcm, drho_pcm);

			double scale_q = (lambda_pcm * m_pcm) / (rho_pcm * dx_pcm*dx_pcm * N3);


			for (svULong j=0; j<np; ++j) {
				dq_dp_adj(i,j) = scale_q * ( (g_adjSens[i])(j,0) - (g_adjSens[i])(j,1) );
			}
		}*/


		// Computation of dq/dp = dq/dT * dT/dp with forward Sensitivities
		Sonic::DMat dq_dp_fwd (nmp, np);
		for (svULong i=0; i<nmp; ++i) {

			double rho_pcm;
			double drho_pcm;

			double T_N1_i = (*g_traj)(i, N1);
			rho_pcm_formula(T_N1_i, rho_pcm, drho_pcm);

			double scale_q = (lambda_pcm * m_pcm) / (rho_pcm * dx_pcm*dx_pcm * N3);

			for (svULong j=0; j<np; ++j) {
				dq_dp_fwd(i,j) = scale_q * ( (g_fwdSens[i])(N1,j) - (g_fwdSens[i])(N1+1,j) );
			}
		}

		// Vergleich fwdSens und adjSens (und speichern in .txt file)
		/*double TOL = 1e-13;
		std::ofstream file_diff_sens;
  		file_diff_sens.open ("/home/argo/masterarbeit/diff_sens.txt");
  		std::ofstream file_sens;
  		file_sens.open ("/home/argo/masterarbeit/sens.txt");
  		file_sens << std::scientific;
  		file_sens.precision(16);
		for (svULong i=0; i<nmp; ++i) {
			for (svULong j=0; j<np; ++j) {
		
				file_sens << (g_fwdSens[i])(N1,j) << "\t" << (g_adjSens[i])(j,0) << "\t" 
					      << (g_fwdSens[i])(N1+1,j) << "\t" << (g_adjSens[i])(j,1) << "\t";

				double relErr1;
				double fwdSens_temp = (g_fwdSens[i])(N1,j);
				double adjSens_temp = (g_adjSens[i])(j,0);
				if ( (fwdSens_temp > TOL && adjSens_temp > TOL) ||
					 (fwdSens_temp > TOL && adjSens_temp < TOL) ||
					 (fwdSens_temp < TOL && adjSens_temp > TOL) ){
					relErr1 = fabs(1 - (g_fwdSens[i])(N1,j) / (g_adjSens[i])(j,0));
				} else if ( fwdSens_temp < TOL || adjSens_temp < TOL ) {
					relErr1 = 0.;
				}
				else {
					relErr1 = 666.666;
				}

				double relErr2;
				fwdSens_temp = (g_fwdSens[i])(N1+1,j);
				adjSens_temp = (g_adjSens[i])(j,1);
				if ( (fwdSens_temp > TOL && adjSens_temp > TOL) ||
					 (fwdSens_temp > TOL && adjSens_temp < TOL) ||
					 (fwdSens_temp < TOL && adjSens_temp > TOL) ){
					relErr2 = fabs(1 - (g_fwdSens[i])(N1,j) / (g_adjSens[i])(j,0));
				} else if ( fwdSens_temp < TOL || adjSens_temp < TOL ) {
					relErr2 = 0.;
				}
				else {
					relErr2 = 666.666;
				}

				file_diff_sens << relErr1 << "\t" << relErr2 << "\t";
			}
			file_diff_sens << std::endl;
			file_sens << std::endl;
		}
		file_diff_sens.close();
		file_sens.close();*/


		// Compute heat flux and residuum vector
		Sonic::DMat heat_flux(nmp, 1);
		for (svULong i=0; i<nmp; ++i) {

			double rho_pcm;
			double drho_pcm;

			double T_N1_i = (*g_traj)(i, N1);
			rho_pcm_formula(T_N1_i, rho_pcm, drho_pcm);

			double scale_q = (lambda_pcm * m_pcm) / (rho_pcm * dx_pcm*dx_pcm * N3);


			heat_flux(i,0) = -scale_q * ((*g_traj)(i,N1+1) - (*g_traj)(i,N1));
		}
		
		Sonic::DMat residuum(nmp, 1);
		for (svULong i=0; i<nmp; ++i) {
			residuum(i,0) = heat_flux(i,0) - (*q_meas)(i,0);
		}



		/******************** Output to Matlab *********************/
		// Residuum Vector
		const Sonic::cDMat& residuum_temp(residuum);
		plhs[0] = mxCreateDoubleMatrix(residuum_temp.nRows(), residuum_temp.nCols(), mxREAL);
		Sonic::DMat residuum_output(mxGetPr(plhs[0]), residuum_temp.nRows(),
				residuum_temp.nRows(), residuum_temp.nCols(), Sonic::CreateReference);
		residuum_output <<= residuum_temp;


		// Jacobian dq/dp
		const Sonic::cDMat& dq_dp_temp(dq_dp_fwd);
		plhs[1] = mxCreateDoubleMatrix(dq_dp_temp.nRows(), dq_dp_temp.nCols(), mxREAL);
		Sonic::DMat dq_dp_output(mxGetPr(plhs[1]), dq_dp_temp.nRows(),
				dq_dp_temp.nRows(), dq_dp_temp.nCols(), Sonic::CreateReference);
		dq_dp_output <<= dq_dp_temp;

		// Jacobian dq/dp
		/*const Sonic::cDMat& dq_dp_temp2(dq_dp_adj);
		plhs[2] = mxCreateDoubleMatrix(dq_dp_temp2.nRows(), dq_dp_temp2.nCols(), mxREAL);
		Sonic::DMat dq_dp_output2(mxGetPr(plhs[2]), dq_dp_temp2.nRows(),
				dq_dp_temp2.nRows(), dq_dp_temp2.nCols(), Sonic::CreateReference);
		dq_dp_output2 <<= dq_dp_temp2;*/


		// Temperature trajectory
		if (nlhs >= 3) {
			const Sonic::cDMat& T_temp(*g_traj);
			plhs[2] = mxCreateDoubleMatrix(T_temp.nRows(), T_temp.nCols(), mxREAL);
			Sonic::DMat T_output(mxGetPr(plhs[2]), T_temp.nRows(),
					T_temp.nRows(), T_temp.nCols(), Sonic::CreateReference);
			T_output <<= T_temp;
		}

		return;


	} else if (strncmp(command, "reset", 99) == 0) {
		// Reset all global variables

		std::cout << "Resetting..." << std::endl;

		delete integrator;
		//delete evaluator;

		nAdjDir = 0;
		nAdjDirTotal = 0;

		delete[] adjSensDir;

		nFwdDir = 0;

		delete[] fwdSensDir;

		solGrid.clear();

		delete q_meas;

		nmp = 0;
		nxd = 0;
		np = 0;

		delete g_traj;
		g_fwdSens.clear();
		g_adjSens.clear();


		// Simulation parameters
		L1 = 0;				
		L3 = 0;				
		N1 = 0;				
		N3 = 0;				
		dx_Const = 0;		
		dx_pcm = 0;			
		alpha = 0;	

		spatial_gridsize.clear();
		alpha_vec.clear();

		lambda_Const = 0;	
		rho_Const = 0;		
		c_p_Const = 0;		
		lambda_pcm = 0;	
		rho_pcm = 0;			

		m_pcm = 0;			

		heat_rate = 0;		
		T_0 = 0;				
		T_end = 0;			

		// Global optimization parameter
		c_p_param_type = uninitialized;
		c_p_params = NULL;

		initialized = false;
	}


}
