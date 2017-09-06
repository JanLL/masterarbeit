#ifndef NURBS_H
#define NURBS_H


#include <iostream>
#include <vector>


class Nurbs
{
	public:
		Nurbs(int num_cntrl_pts, int nurbs_order);

		std::vector<double> get_cntrl_pts_x();
		std::vector<double> get_cntrl_pts_y();
		std::vector<double> get_weights();
		std::vector<double> get_U();
		
		void set_cntrl_pts_x(std::vector<double> cntrl_pts_x);
		void set_cntrl_pts_y(std::vector<double> cntrl_pts_y);
		//void set_U(std::vector<double> U);
		//void set_weights(std::vector<double> weights);

		double eval_basis_fcn(double u, int i, int p);


		std::vector<double> eval_nurbs_curve(double u);

		int get_interval_index(double u);

		double compute_a_i_m3(int i, double u);
		double compute_a_i_m2(int i, double u);
		double compute_a_i_m1(int i, double u);
		double compute_a_i_0(int i, double u);

		std::vector<double> eval_nurbs_curve_mod(double u);



	private:
		int num_cntrl_pts;
		int nurbs_order;
		int len_U;

		std::vector<double> cntrl_pts_x;
		std::vector<double> cntrl_pts_y;
		std::vector<double> weights;
		std::vector<double> U;

				
};



#endif