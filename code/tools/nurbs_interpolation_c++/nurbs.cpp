#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>

#include "interpolation.h"

#include "nurbs.hpp"

#include <adolc/adouble.h>


template<typename T>
Nurbs<T>::Nurbs(int num_cntrl_pts_in, int nurbs_order_in) :
	num_cntrl_pts(num_cntrl_pts_in), 
	nurbs_order(nurbs_order_in),
	len_U(num_cntrl_pts + nurbs_order),
	cntrl_pts_x(num_cntrl_pts_in, 0),
	cntrl_pts_y(num_cntrl_pts_in, 0),
	weights(num_cntrl_pts_in, 1.),
	U(num_cntrl_pts_in + nurbs_order_in)
{

	for (int i=0; i <= nurbs_order; ++i) {
		U[i] = 0.;
		U[len_U-1-i] = 1.;
	}

	for (int i=1; i <= num_cntrl_pts - nurbs_order; ++i) {
		U[i + nurbs_order - 1] = i / double(num_cntrl_pts - nurbs_order + 1);
	}
}


template<typename T>
std::vector<T> Nurbs<T>::get_cntrl_pts_x() {
	return cntrl_pts_x;
}

template<typename T>
std::vector<T> Nurbs<T>::get_cntrl_pts_y() {
	return cntrl_pts_y;
}

template<typename T>
std::vector<T> Nurbs<T>::get_weights() {
	return weights;
}

template<typename T>
std::vector<T> Nurbs<T>::get_U() {
	return U;
}

template<typename T>
void Nurbs<T>::set_cntrl_pts_x(std::vector<T> cntrl_pts_x_in) {
	for (int i=0; i <= num_cntrl_pts-1; ++i) {
		cntrl_pts_x[i] = cntrl_pts_x_in[i];
	}
}
template<typename T>
void Nurbs<T>::set_cntrl_pts_y(std::vector<T> cntrl_pts_y_in) {
	for (int i=0; i <= num_cntrl_pts-1; ++i) {
		cntrl_pts_y[i] = cntrl_pts_y_in[i];
	}
}


template<typename T>
T Nurbs<T>::eval_basis_fcn(T u, int i, int p) {
	if (p == 0) {
		if (u >= U[i] && u <= U[i+1]) {
			return 1.;
		}
		else {
			return 0.;
		}
	}
	else {

		T part1 = 0;
		T part2 = 0;

		// if denominator is zero, N is also zero, set then part1=0
		T N_i_pm1 = this->eval_basis_fcn(u, i, p-1);
		if (N_i_pm1 != 0) {
			part1 = (u - U[i]) / (U[i+p] - U[i]) * N_i_pm1;
		}

		T N_ip1_pm1 = this->eval_basis_fcn(u, i+1, p-1);
		if (N_ip1_pm1 != 0) {
			part2 = (U[i+p+1] - u) / (U[i+p+1] - U[i+1]) * N_ip1_pm1;
		}

		return part1 + part2;
	}
}


template<typename T>
std::vector<T> Nurbs<T>::eval_nurbs_curve(T u) {

	std::vector<T> N_vec(num_cntrl_pts);


	for (int i=0; i < num_cntrl_pts; ++i) {
		N_vec[i] = this->eval_basis_fcn(u, i, nurbs_order-1);
	}

	
	T denominator = 0;
	T numerator_x = 0;
	T numerator_y = 0;

	for (int i=0; i < num_cntrl_pts; ++i) {
		numerator_x += N_vec[i] * weights[i] * cntrl_pts_x[i];
		numerator_y += N_vec[i] * weights[i] * cntrl_pts_y[i];
		denominator += N_vec[i] * weights[i];

	}

	std::vector<T> C(2);
	C[0] = numerator_x / denominator;
	C[1] = numerator_y / denominator;

	return C;
}


template<typename T>
T Nurbs<T>::compute_a_i_m3(int i, T u) {

	T a_i_m3;

	a_i_m3 = (U[i+4] - u)*(U[i+4] - u)*(U[i+4] - u) / ((U[i+4] - U[i+1])*(U[i+4] - U[i+2])*(U[i+4] - U[i+3]));

	return a_i_m3;
}


template<typename T>
T Nurbs<T>::compute_a_i_m2(int i, T u) {

	T a_i_m2;

	a_i_m2 = (u - U[i])/(U[i+3] - U[i]) * (U[i+3] - u)/(U[i+3] - U[i+1]) * (U[i+3] - u)/(U[i+3] - U[i+2]) + 
			 (U[i+4] - u)/(U[i+4] - U[i+1]) * ( (u - U[i+1])/(U[i+3] - U[i+1]) * (U[i+3] - u)/(U[i+3] - U[i+2])
			                                   +(U[i+4] - u)/(U[i+4] - U[i+2]) * (u - U[i+2])/(U[i+3] - U[i+2]) );

	return a_i_m2;
}


template<typename T>
T Nurbs<T>::compute_a_i_m1(int i, T u) {

	T a_i_m1;

	a_i_m1 = (u - U[i])/(U[i+3] - U[i]) * ( (u - U[i])/(U[i+2] - U[i]) * (U[i+2] - u)/(U[i+2] - U[i+1])
			                               +(U[i+3] - u)/(U[i+3] - U[i+1]) * (u - U[i+1])/(U[i+2] - U[i+1]) )
	         + (U[i+4] - u)/(U[i+4] - U[i+1]) * (u - U[i+1])/(U[i+3] - U[i+1]) * (u - U[i+1])/(U[i+2] - U[i+1]);

	return a_i_m1;
}


template<typename T>
T Nurbs<T>::compute_a_i_0(int i, T u) {

	T a_i_0;

	a_i_0 = (u - U[i])*(u - U[i])*(u - U[i]) / ((U[i+3] - U[i])*(U[i+2] - U[i])*(U[i+1] - U[i]));

	return a_i_0;
}


template<typename T>
int Nurbs<T>::get_interval_index(T u) {

	int i=3;  // NURBS order = 4
	for (int j=3; j <= len_U-4; ++j) {

		if (u - U[j] > 0) {
			i = j;
		}
		else {
			break;
		}

	}

	//T du = U[nurbs_order] - U[nurbs_order-1];
	//int i = floor(u/du) + 3;

	return i;
}


template<typename T>
T Nurbs<T>::eval_nurbs_curve_x(T u) {

	int i = this->get_interval_index(u);

	T a_im3_m3 = this->compute_a_i_m3(i-3, u);
	T a_im2_m2 = this->compute_a_i_m2(i-2, u);
	T a_im1_m1 = this->compute_a_i_m1(i-1, u);
	T a_i_0    = this->compute_a_i_0 (i  , u);


	T numerator_x = a_im3_m3 * weights[i-3] * cntrl_pts_x[i-3] +
						 a_im2_m2 * weights[i-2] * cntrl_pts_x[i-2] +
						 a_im1_m1 * weights[i-1] * cntrl_pts_x[i-1] +
						 a_i_0    * weights[i-0] * cntrl_pts_x[i-0];

	T denominator = a_im3_m3 * weights[i-3] +
						 a_im2_m2 * weights[i-2] +
						 a_im1_m1 * weights[i-1] +
						 a_i_0    * weights[i-0];

	return numerator_x / denominator;
}	


template<typename T>
T Nurbs<T>::eval_nurbs_curve_y(T u) {

	int i = this->get_interval_index(u);

	T a_im3_m3 = this->compute_a_i_m3(i-3, u);
	T a_im2_m2 = this->compute_a_i_m2(i-2, u);
	T a_im1_m1 = this->compute_a_i_m1(i-1, u);
	T a_i_0 = this->compute_a_i_0(i, u);

	T numerator_y = a_im3_m3 * weights[i-3] * cntrl_pts_y[i-3] +
						 a_im2_m2 * weights[i-2] * cntrl_pts_y[i-2] +
						 a_im1_m1 * weights[i-1] * cntrl_pts_y[i-1] +
						 a_i_0    * weights[i-0] * cntrl_pts_y[i-0];

	T denominator = a_im3_m3 * weights[i-3] +
						 a_im2_m2 * weights[i-2] +
						 a_im1_m1 * weights[i-1] +
						 a_i_0    * weights[i-0];

	return numerator_y / denominator;
}


template<typename T>
T Nurbs<T>::get_u_from_Cx(T Cx, T u_start, T TOL) {

	T err = 1e8;
	T u = u_start;
	T Cx_newton;
	T du;
	T h = 0.005; 

	while (err > TOL) {

		Cx_newton = this->eval_nurbs_curve_x(u);

		du = -(Cx_newton - Cx);

		u += h*du;

		err = abs(Cx_newton - Cx);
	}

	return u;
}




int main(int argc, char** argv) {


	int nurbs_order = 4;
	int poly_order = nurbs_order - 1;

	std::vector<adouble> cntrl_pts_x = {0, 30, 60, 90, 120, 125, 130, 132., 135, 150, 160, 180};
	std::vector<adouble> cntrl_pts_y = {1., 1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52, 1.53, 1.54};
	

	int num_cntrl_pts = cntrl_pts_x.size();


	Nurbs<adouble> nurbs(num_cntrl_pts, nurbs_order);
	nurbs.set_cntrl_pts_x(cntrl_pts_x);
	nurbs.set_cntrl_pts_y(cntrl_pts_y);


	std::ofstream file_nurbes;
  	file_nurbes.open ("curve_nurbes.txt");

  	for (adouble u=0.; u<1; u+=0.001) {

		adouble Cx = nurbs.eval_nurbs_curve_x(u);
		adouble Cy = nurbs.eval_nurbs_curve_y(u);
		
		file_nurbes << u << "\t" << Cx << "\t" << Cy << std::endl;

	}
	file_nurbes.close();

	return 0;
	}/*

	// TODO: Newton implementieren und Laufzeittest
	double TOL = 0.1;
	double u_start = 0.;

  	double u = u_start;

  	for (double Ts=50.; Ts < 55.; Ts+=0.1) {
		clock_t begin = clock();
		
		u = nurbs.get_u_from_Cx(Ts, u, TOL);

		clock_t end = clock();
		std::cout << std::setprecision(20) << double(end - begin) / CLOCKS_PER_SEC << std::endl;
	}
	//clock_t end = clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//std::cout << std::setprecision(20) << elapsed_secs << std::endl;


	}
	

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