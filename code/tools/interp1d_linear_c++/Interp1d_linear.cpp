#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "Interp1d_linear.hpp"




Interp1d_linear::Interp1d_linear(std::vector<double> x_data_input, std::vector<double> y_data_input) {

	for(std::vector<double>::iterator it=x_data_input.begin(); it!=x_data_input.end(); ++it) {
		x_data.insert(*it);
	}

	for(std::vector<double>::iterator it=y_data_input.begin(); it!=y_data_input.end(); ++it) {
		y_data.insert(*it);
	}

	compute_coeffs();

	
}

void Interp1d_linear::compute_coeffs() {

	double x_i, x_ip1, y_i, y_ip1;
	std::vector<double> coeffs_i(2);
	std::set<double>::iterator it_x_next;
	std::set<double>::iterator it_y_next;

	for (std::set<double>::iterator it_x=x_data.begin(), it_y=y_data.begin(); it_x!=x_data.end(); ++it_x, ++it_y) {

		// Nicht ganz optimal, iterator auf naechstes element koennte man in naechster iteration der for-loop wieder verwenden...
		it_x_next = it_x;
		it_x_next++;
		it_y_next = it_y;
		it_y_next++;

		x_i = *it_x;
		x_ip1 = *it_x_next;
		y_i = *it_y;
		y_ip1 = *it_y_next;

		coeffs_i[1] = (y_ip1 - y_i) / (x_ip1 - x_i);
		coeffs_i[0] = y_i - coeffs_i[1]*x_i;

		coeffs_map[*it_x] = coeffs_i;

	}	


}

double Interp1d_linear::operator()(double x) {

	double y_interp;
	std::set<double>::iterator it_i;
	std::vector<double> coeffs(2);

	it_i = x_data.lower_bound(x);
	it_i--;

	coeffs = coeffs_map[*it_i];

	y_interp = coeffs[0] + coeffs[1]*x;

	return y_interp;

}





int main(int argc, char** argv) {

	std::cout << "blubb!\n";

	std::vector<double> x_data(5);
	std::vector<double> y_data(5);

	for(int i=0; i<5; ++i) {
		x_data[i] = i+1;
		y_data[i] = i*i;
	}

	Interp1d_linear interpolator(x_data, y_data);

	/*std::vector<double> coeff_pair(2, 1);
	for(int i=0; i<5; ++i) {
		coeff_pair = interpolator.coeffs_map[i];
		std::cout << coeff_pair[0] << "\t" << coeff_pair[1] << std::endl;
	}*/

	std::cout << interpolator(3.5) << std::endl;




	return 0;
}