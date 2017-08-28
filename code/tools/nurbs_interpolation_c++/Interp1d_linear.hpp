#ifndef INTERP1D_LINEAR_H
#define INTERP1D_LINEAR_H


#include <iostream>
#include <vector>
#include <set>
#include <map>



template<typename T>
class Interp1d_linear {

	public:
		Interp1d_linear(std::vector<T> x_data, std::vector<T> y_data);

		std::vector<T> operator()(T x);

		T eval(T x);
		T eval_d(T x);

	private:
		void compute_coeffs();
		std::vector<T> get_coeffs(T x);


		std::vector<T> x_data;
		std::vector<T> y_data;

		int num_data_pts;


		std::map<T, std::vector<T> > coeffs_map;


};


#endif