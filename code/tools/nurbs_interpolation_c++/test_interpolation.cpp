#include <iostream>
#include <fstream>

#include "Interp1d_linear.hpp"




int main(int argc, char** argv) {

	std::vector<double> x_data(10);
	std::vector<double> y_data(10);

	for(int i=0; i<10; ++i) {
		x_data[i] = i;
	}

	y_data[0] = 1.;
	y_data[1] = 4.;
	y_data[2] = 3.;
	y_data[3] = 7.;
	y_data[4] = 2.;
	y_data[5] = 2.;
	y_data[6] = -1.;
	y_data[7] = 3.;
	y_data[8] = 4.;
	y_data[9] = 2.;

	Interp1d_linear interpolator(x_data, y_data);


	std::ofstream file;
	file.open("data_interp.txt");

	double h = 0.01;
	for (double d=0.+h; d <= 9.; d += h) {
		double y_interp = interpolator(d);
		file << d << "\t" << y_interp << "\n";
	}

	return 0;
}