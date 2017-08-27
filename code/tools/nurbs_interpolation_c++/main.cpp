#include <iostream>
#include <fstream>

#include <time.h>


#include "nurbs.hpp"
#include "Interp1d_linear.hpp"


int main(int argc, char** argv) {

	// NURBS part
	int nurbs_order = 4;
	int poly_order = nurbs_order - 1;

	std::vector<double> cntrl_pts_x(12);// = {0, 30, 60, 90, 120, 125, 130, 132., 135, 150, 160, 180};
	std::vector<double> cntrl_pts_y(12);// = {1., 1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52, 1.53, 1.54};
	int num_cntrl_pts = cntrl_pts_x.size();

	cntrl_pts_x[0] = 0; 	cntrl_pts_y[0] = 1.;
	cntrl_pts_x[1] = 30; 	cntrl_pts_y[1] = 1.;
	cntrl_pts_x[2] = 60; 	cntrl_pts_y[2] = 1.1;
	cntrl_pts_x[3] = 90; 	cntrl_pts_y[3] = 1.15;
	cntrl_pts_x[4] = 120; 	cntrl_pts_y[4] = 1.2;
	cntrl_pts_x[5] = 125; 	cntrl_pts_y[5] = 5.;
	cntrl_pts_x[6] = 130; 	cntrl_pts_y[6] = 10.;
	cntrl_pts_x[7] = 132; 	cntrl_pts_y[7] = 1.5;
	cntrl_pts_x[8] = 135; 	cntrl_pts_y[8] = 1.51;
	cntrl_pts_x[9] = 150; 	cntrl_pts_y[9] = 1.52;
	cntrl_pts_x[10] = 160; 	cntrl_pts_y[10] = 1.53;
	cntrl_pts_x[11] = 180; 	cntrl_pts_y[11] = 1.54;
	


	Nurbs nurbs(num_cntrl_pts, nurbs_order);
	nurbs.set_cntrl_pts_x(cntrl_pts_x);
	nurbs.set_cntrl_pts_y(cntrl_pts_y);


	double h = 0.001;

	
	std::vector<double> C_x(int(1/h)+1);
	std::vector<double> C_y(int(1/h)+1);
	

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


	// 1D Interpolation part
	Interp1d_linear interpolator(C_x, C_y);

	interpolator(30.);



	std::ofstream file_interp;
	file_interp.open("curve_interp.txt");

	double time_start, time_duration;

	time_start = clock();  // start time measurement
	for(double d=30.; d <= 160; d+=0.01) {
		std::vector<double> values = interpolator(d);
		file_interp << d << "\t" << values[0] << "\t" << values[1] << std::endl;
	}
	time_duration = (clock() - time_start) / CLOCKS_PER_SEC;

	file_interp.close();

	std::cout << time_duration << " s\n";

	return 0;
}