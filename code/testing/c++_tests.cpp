#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <cmath>

#include <time.h>

#include <adolc/adouble.h>



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

	std::ofstream rho_test;
  	rho_test.open ("/home/argo/masterarbeit/rho_test.txt");

  	double rho;
  	double drho;

  	for (double d=30.; d<200.; d += 0.01) {
  		rho_pcm_formula(d, rho, drho);
  		rho_test << d << "\t" << rho << "\t" << drho << std::endl;
  	}

  	rho_test.close();


	return 0;

}