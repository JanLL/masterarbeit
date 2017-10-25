#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <cmath>

#include <time.h>

#include <adolc/adouble.h>

#include "/home/argo/masterarbeit/code/c_p_parametrizations.cpp"



template<typename T>
T fraser_suzuki(T x, T h, T r, T wr, T sr, T z) {

	T value;

	T log_value = log(1 + (x-z)*(sr*sr-1)/(wr*sr));

	return h*exp(-log(r)/(log(sr)*log(sr)) * log_value*log_value);

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






int main(int argc, char** argv) {

	/*std::ofstream rho_test;
  	rho_test.open ("/home/argo/masterarbeit/rho_test.txt");

  	double rho;
  	double drho;

  	for (double d=30.; d<200.; d += 0.01) {
  		rho_pcm_formula(d, rho, drho);
  		rho_test << d << "\t" << rho << "\t" << drho << std::endl;
  	}

  	rho_test.close();*/


	double p_optim_end[32] = {16.6062086431708,	1.92927098591158,	128.850918977868,	
							  34.1327529190694,	1.27007399531043,	130.601505840094,	
							  4.66622071658818,	3.08315397789580,	126.387021606664,	
							  0.330799303747586,	6.28205711993873,	138.931914796423,	
							 -1.99975880363326,	1.00794170224404,	127.886002738954,	
							  2.00575452880237,	5.23818264316667,	122.475803758605,	
							  8.69795539050226,	0.764167761467260,	131.898468412807,
							  0.,				26.4978102269972,	286.560297948564,	
							  0.,				11.4117833489350,	111.924264531493,
							  0.,				35.7646524748368,	95.2324216508870,	
						      0.,				0.};


	std::ofstream enthalpy_test;
  	enthalpy_test.open ("/home/argo/masterarbeit/enthalpy_test.txt");

	double c_p;
	double h;
	
	for (double T=30; T<200; T+=0.05) {
	
		gauss_linear_comb_formula(T, c_p, h, p_optim_end);

		enthalpy_test << T << "\t" << c_p << "\t" << h << std::endl;

	}

	enthalpy_test.close();

	return 0;

}