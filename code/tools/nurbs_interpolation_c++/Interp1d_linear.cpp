#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <list>

#include "Interp1d_linear.hpp"



template<typename T>
Interp1d_linear<T>::Interp1d_linear(std::vector<T> x_data_input, std::vector<T> y_data_input) : 
	x_data(x_data_input), 
	y_data(y_data_input), 
	num_data_pts(x_data_input.size()) 
{
	compute_coeffs();
}

template<typename T>
void Interp1d_linear<T>::compute_coeffs() {

	T x_i, x_ip1, y_i, y_ip1;
	std::vector<T> coeffs_i(2);

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

template<typename T>
std::vector<T> Interp1d_linear<T>::get_coeffs(T x) {

	typename std::map<T, std::vector<T> >::iterator it_i;

	it_i = coeffs_map.lower_bound(x);
	it_i--;

	std::vector<T> coeffs = (*it_i).second;

	return coeffs;
}


template<typename T>
T Interp1d_linear<T>::eval(T x) {

	T y_interp;

	std::vector<T> coeffs = get_coeffs(x);

	y_interp = coeffs[0] + coeffs[1]*x;

	return y_interp;
}


template<typename T>
T Interp1d_linear<T>::eval_d(T x) {

	T yd_interp;

	std::vector<T> coeffs = get_coeffs(x);

	yd_interp = coeffs[1];

	return yd_interp;
}


template<typename T>
std::vector<T> Interp1d_linear<T>::operator()(T x) {

	// Note: Currently x must be: x_data[0] < x <= x_data[end]


	std::vector<T> values(2);

	std::vector<T> coeffs = get_coeffs(x);


	values[0] = coeffs[0] + coeffs[1]*x;
	values[1] = coeffs[1];


	return values;

}





