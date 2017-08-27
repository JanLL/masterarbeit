#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <list>

#include "Interp1d_linear.hpp"




Interp1d_linear::Interp1d_linear(std::vector<double> x_data_input, std::vector<double> y_data_input) : 
	x_data(x_data_input), 
	y_data(y_data_input), 
	num_data_pts(x_data_input.size()) 
{
	compute_coeffs();
}

void Interp1d_linear::compute_coeffs() {

	double x_i, x_ip1, y_i, y_ip1;
	std::vector<double> coeffs_i(2);

	for (int i=0; i < num_data_pts-1; ++i) {

		x_i   = x_data[i];
		x_ip1 = x_data[i+1];
		y_i   = y_data[i];
		y_ip1 = y_data[i+1];

		coeffs_i[1] = (y_ip1 - y_i) / (x_ip1 - x_i); // slope
		coeffs_i[0] = y_i - coeffs_i[1]*x_i;         // offset

		coeffs_map[x_i] = coeffs_i;


	}	

}


std::vector<double> Interp1d_linear::get_coeffs(double x) {

	std::map<double, std::vector<double> >::iterator it_i;

	it_i = coeffs_map.lower_bound(x);
	it_i--;

	std::vector<double> coeffs = (*it_i).second;

	return coeffs;
}



double Interp1d_linear::eval(double x) {

	double y_interp;

	std::vector<double> coeffs = get_coeffs(x);

	y_interp = coeffs[0] + coeffs[1]*x;

	return y_interp;
}


double Interp1d_linear::eval_d(double x) {

	double yd_interp;

	std::vector<double> coeffs = get_coeffs(x);

	yd_interp = coeffs[1];

	return yd_interp;
}


std::vector<double> Interp1d_linear::operator()(double x) {

	// Note: Currently x must be: x_data[0] < x <= x_data[end]


	std::vector<double> values(2);

	std::vector<double> coeffs = get_coeffs(x);


	values[0] = coeffs[0] + coeffs[1]*x;
	values[1] = coeffs[1];


	return values;

}





