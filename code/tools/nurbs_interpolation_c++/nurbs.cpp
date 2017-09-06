#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>

#include "interpolation.h"

#include "nurbs.hpp"

//#include <adolc/adouble.h>



Nurbs::Nurbs(int num_cntrl_pts_in, int nurbs_order_in) :
	num_cntrl_pts(num_cntrl_pts_in), 
	nurbs_order(nurbs_order_in),
	len_U(num_cntrl_pts + nurbs_order),
	cntrl_pts_x(num_cntrl_pts_in, 0),
	cntrl_pts_y(num_cntrl_pts_in, 0),
	weights(num_cntrl_pts_in, 1.),
	U(num_cntrl_pts_in + nurbs_order_in)
{

	for (int i=0; i <= nurbs_order; ++i) {
		U[i] = 0;
		U[len_U-1-i] = 1;
	}

	for (int i=1; i <= num_cntrl_pts - nurbs_order; ++i) {
		U[i + nurbs_order - 1] = i / float(num_cntrl_pts - nurbs_order + 1);
	}
}




std::vector<double> Nurbs::get_cntrl_pts_x() {
	return cntrl_pts_x;
}

std::vector<double> Nurbs::get_cntrl_pts_y() {
	return cntrl_pts_y;
}

std::vector<double> Nurbs::get_weights() {
	return weights;
}

std::vector<double> Nurbs::get_U() {
	return U;
}


void Nurbs::set_cntrl_pts_x(std::vector<double> cntrl_pts_x_in) {
	for (int i=0; i <= num_cntrl_pts-1; ++i) {
		cntrl_pts_x[i] = cntrl_pts_x_in[i];
	}
}

void Nurbs::set_cntrl_pts_y(std::vector<double> cntrl_pts_y_in) {
	for (int i=0; i <= num_cntrl_pts-1; ++i) {
		cntrl_pts_y[i] = cntrl_pts_y_in[i];
	}
}



double Nurbs::eval_basis_fcn(double u, int i, int p) {
	if (p == 0) {
		if (u >= U[i] && u <= U[i+1]) {
			return 1.;
		}
		else {
			return 0.;
		}
	}
	else {

		double part1 = 0;
		double part2 = 0;

		// if denominator is zero, N is also zero, set then part1=0
		double N_i_pm1 = this->eval_basis_fcn(u, i, p-1);
		if (N_i_pm1 != 0) {
			part1 = (u - U[i]) / (U[i+p] - U[i]) * N_i_pm1;
		}

		double N_ip1_pm1 = this->eval_basis_fcn(u, i+1, p-1);
		if (N_ip1_pm1 != 0) {
			part2 = (U[i+p+1] - u) / (U[i+p+1] - U[i+1]) * N_ip1_pm1;
		}

		return part1 + part2;
	}
}


std::vector<double> Nurbs::eval_nurbs_curve(double u) {

	std::vector<double> N_vec(num_cntrl_pts);


	for (int i=0; i < num_cntrl_pts; ++i) {
		N_vec[i] = this->eval_basis_fcn(u, i, nurbs_order-1);
	}

	
	double denominator = 0;
	double numerator_x = 0;
	double numerator_y = 0;

	for (int i=0; i < num_cntrl_pts; ++i) {
		numerator_x += N_vec[i] * weights[i] * cntrl_pts_x[i];
		numerator_y += N_vec[i] * weights[i] * cntrl_pts_y[i];
		denominator += N_vec[i] * weights[i];

	}

	std::vector<double> C(2);
	C[0] = numerator_x / denominator;
	C[1] = numerator_y / denominator;

	return C;

}

double Nurbs::compute_a_i_m3(int i, double u) {

	double a_i_m3;

	a_i_m3 = (U[i+4] - u)*(U[i+4] - u)*(U[i+4] - u) / ((U[i+4] - U[i+1])*(U[i+4] - U[i+2])*(U[i+4] - U[i+3]));

	return a_i_m3;

}


double Nurbs::compute_a_i_m2(int i, double u) {

	double a_i_m2;

	a_i_m2 = (u - U[i])/(U[i+3] - U[i]) * (U[i+3] - u)/(U[i+3] - U[i+1]) * (U[i+3] - u)/(U[i+3] - U[i+2]) + 
			 (U[i+4] - u)/(U[i+4] - U[i+1]) * ( (u - U[i+1])/(U[i+3] - U[i+1]) * (U[i+3] - u)/(U[i+3] - U[i+2])
			                                   +(U[i+4] - u)/(U[i+4] - U[i+2]) * (u - U[i+2])/(U[i+3] - U[i+2]) );

	return a_i_m2;

}

double Nurbs::compute_a_i_m1(int i, double u) {

	double a_i_m1;

	a_i_m1 = (u - U[i])/(U[i+3] - U[i]) * ( (u - U[i])/(U[i+2] - U[i]) * (U[i+2] - u)/(U[i+2] - U[i+1])
			                               +(U[i+3] - u)/(U[i+3] - U[i+1]) * (u - U[i+1])/(U[i+2] - U[i+1]) )
	         + (U[i+4] - u)/(U[i+4] - U[i+1]) * (u - U[i+1])/(U[i+3] - U[i+1]) * (u - U[i+1])/(U[i+2] - U[i+1]);

	return a_i_m1;

}

double Nurbs::compute_a_i_0(int i, double u) {

	double a_i_0;

	a_i_0 = (u - U[i])*(u - U[i])*(u - U[i]) / ((U[i+3] - U[i])*(U[i+2] - U[i])*(U[i+1] - U[i]));

	return a_i_0;

}


int Nurbs::get_interval_index(double u) {

	double du = U[nurbs_order] - U[nurbs_order-1];
	int i = floor(u/du) + 3;

	return i;

}


std::vector<double> Nurbs::eval_nurbs_curve_mod(double u) {

	int i = this->get_interval_index(u);

	double a_im3_m3 = this->compute_a_i_m3(i-3, u);
	double a_im2_m2 = this->compute_a_i_m2(i-2, u);
	double a_im1_m1 = this->compute_a_i_m1(i-1, u);
	double a_i_0 = this->compute_a_i_0(i, u);


	double numerator_x = a_im3_m3 * weights[i-3] * cntrl_pts_x[i-3] +
						 a_im2_m2 * weights[i-2] * cntrl_pts_x[i-2] +
						 a_im1_m1 * weights[i-1] * cntrl_pts_x[i-1] +
						 a_i_0    * weights[i-0] * cntrl_pts_x[i-0];

	double numerator_y = a_im3_m3 * weights[i-3] * cntrl_pts_y[i-3] +
						 a_im2_m2 * weights[i-2] * cntrl_pts_y[i-2] +
						 a_im1_m1 * weights[i-1] * cntrl_pts_y[i-1] +
						 a_i_0    * weights[i-0] * cntrl_pts_y[i-0];

	double denominator = a_im3_m3 * weights[i-3] +
						 a_im2_m2 * weights[i-2] +
						 a_im1_m1 * weights[i-1] +
						 a_i_0    * weights[i-0];

	std::vector<double> C(2);
	C[0] = numerator_x / denominator;
	C[1] = numerator_y / denominator;

	return C;

}	




int main(int argc, char** argv) {

	int nurbs_order = 4;
	int poly_order = nurbs_order - 1;

	std::vector<double> cntrl_pts_x = {0, 30, 60, 90, 120, 125, 130, 132., 135, 150, 160, 180};
	std::vector<double> cntrl_pts_y = {1., 1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52, 1.53, 1.54};
	//std::vector<double> cntrl_pts_x = {0, 1, 2, 3, 4, 5, 6, 7};
	//std::vector<double> cntrl_pts_y = {0, 1, 2, 3, 4, 5, 6, 7};
	

	int num_cntrl_pts = cntrl_pts_x.size();


	Nurbs nurbs(num_cntrl_pts, nurbs_order);
	nurbs.set_cntrl_pts_x(cntrl_pts_x);
	nurbs.set_cntrl_pts_y(cntrl_pts_y);

	/*std::vector<double> U = nurbs.get_U();
	for (std::vector<double>::iterator it=U.begin(); it!=U.end(); ++it) {
		std::cout << *it << "\t";
	}
	std::cout << std::endl;*/


	std::ofstream file_nurbes;
  	file_nurbes.open ("curve_nurbes.txt");

  	for (double u=0.; u<1; u+=0.001) {

		std::vector<double> C = nurbs.eval_nurbs_curve(u);
		std::vector<double> D = nurbs.eval_nurbs_curve_mod(u);
		
		file_nurbes << C[0] << "\t" << C[1] << "\t" << D[0] << "\t" << D[1] << std::endl;

	}
	file_nurbes.close();

	}/*
	

	double h = 0.01;

	alglib::real_1d_array C_x;
	C_x.setlength(int(1/h) + 1);
	alglib::real_1d_array C_y;
	C_y.setlength(int(1/h) + 1);
	

	std::ofstream file_nurbes;
  	file_nurbes.open ("curve_nurbes.txt");

	std::vector<double> C_temp(2);
	int i=0;
	for (double u=0. + 1e-8; u <= 1.; u+=h) {
		C_temp = nurbs.eval_nurbs_curve(u);
		C_x[i] = C_temp[0];
		C_y[i] = C_temp[1];

		file_nurbes << C_x[i] << "\t" << C_y[i] << std::endl;

		i++;
	}
	file_nurbes.close();

	

	alglib::spline1dinterpolant interpolant; 
	alglib::spline1dbuildlinear(C_x, C_y, interpolant);



	double c_p;
	double dc_p;
	double ddc_p;

	std::ofstream file_interpol;
  	file_interpol.open ("curve_interpol.txt");

	double time_start, time_duration;

	time_start = clock();  // start time measurement
	for (double T=30.; T <= 160; T += 0.01) {	

		alglib::spline1ddiff(interpolant, T, c_p, dc_p, ddc_p);
		//file << T << "\t" << c_p << "\t" << dc_p << "\t" << ddc_p << std::endl; 
		file_interpol << T << "\t" << alglib::spline1dcalc(interpolant, T) << "\n";
	}
	time_duration = (clock() - time_start) / CLOCKS_PER_SEC;
	file_interpol.close();

	std::cout << time_duration << std::endl;

	return 0;

}*/