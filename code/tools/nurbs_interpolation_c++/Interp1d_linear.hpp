#include <iostream>
#include <vector>
#include <set>
#include <map>




class Interp1d_linear {

	public:
		Interp1d_linear(std::vector<double> x_data, std::vector<double> y_data);

		std::vector<double> operator()(double x);

		double eval(double x);
		double eval_d(double x);

	private:
		void compute_coeffs();
		std::vector<double> get_coeffs(double x);

		int num_data_pts;

		std::vector<double> x_data;
		std::vector<double> y_data;

		std::map<double, std::vector<double> > coeffs_map;


};