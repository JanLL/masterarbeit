#include <iostream>
#include <vector>
#include <set>
#include <map>




class Interp1d_linear {

	public:
		Interp1d_linear(std::vector<double> x_data, std::vector<double> y_data);

		double operator()(double x);

	//private:
		void compute_coeffs();

		std::set<double> x_data;
		std::set<double> y_data;

		std::map<double, std::vector<double> > coeffs_map;


};