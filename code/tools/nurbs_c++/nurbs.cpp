#include <iostream>
#include <vector>


#include "nurbs.hpp"


Nurbs::Nurbs(int num_cntrl_pts_in, int nurbs_order_in) :
	num_cntrl_pts(num_cntrl_pts_in), 
	nurbs_order(nurbs_order_in),
	cntrl_pts_x(num_cntrl_pts_in, 0),
	cntrl_pts_y(num_cntrl_pts_in, 0),
	weights(num_cntrl_pts_in, 1.),
	U(num_cntrl_pts_in + nurbs_order_in)
{
	int len_U = num_cntrl_pts + nurbs_order;

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
	//std::cout << "eval_basis_fcn called\n";
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


// TODO: Funktionen zum Auswerten der NURBES Kurve ...



int main(int argc, char** argv) {

	int num_cntrl_pts = 10;
	int nurbs_order = 4;

	int poly_order = nurbs_order - 1;

	Nurbs nurbs(num_cntrl_pts, nurbs_order);

	std::vector<double> cntrl_pts_x = {1., 2., 3., 4., 5., 6., 7., 8.};
	std::vector<double> cntrl_pts_y = {2., 4., 1., 7., 10., 2., 1., 3.};
	
	nurbs.set_cntrl_pts_x(cntrl_pts_x);
	nurbs.set_cntrl_pts_y(cntrl_pts_y);

	double u = 0.5;
	
	for (int i=0; i <= num_cntrl_pts; ++i) {
		std::cout << nurbs.eval_basis_fcn(0.5, i, nurbs_order) << std::endl;
	}


	return 0;

}