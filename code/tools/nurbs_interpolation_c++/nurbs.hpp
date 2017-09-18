#ifndef NURBS_H
#define NURBS_H


#include <iostream>
#include <vector>

template<typename T>
class Nurbs
{
	public:
		Nurbs(int num_cntrl_pts, int nurbs_order);

		std::vector<T> get_cntrl_pts_x();
		std::vector<T> get_cntrl_pts_y();
		std::vector<T> get_weights();
		std::vector<T> get_U();
		
		void set_cntrl_pts_x(std::vector<T> cntrl_pts_x);
		void set_cntrl_pts_y(std::vector<T> cntrl_pts_y);
		//void set_U(std::vector<T> U);
		//void set_weights(std::vector<T> weights);

		T eval_basis_fcn(T u, int i, int p);


		std::vector<T> eval_nurbs_curve(T u);

		int get_interval_index(T u);

		T compute_a_i_m3(int i, T u);
		T compute_a_i_m2(int i, T u);
		T compute_a_i_m1(int i, T u);
		T compute_a_i_0(int i, T u);


		T compute_da_i_m3(int i, T u);
		T compute_da_i_m2(int i, T u);
		T compute_da_i_m1(int i, T u);
		T compute_da_i_0(int i, T u);

		T eval_nurbs_curve_x(T u);
		T eval_nurbs_curve_y(T u);

		T eval_dCy_dCx(T u);

		T get_u_from_Cx(T Cx, T u_start, T TOL);




	private:
		int num_cntrl_pts;
		int nurbs_order;
		int len_U;

		std::vector<T> cntrl_pts_x;
		std::vector<T> cntrl_pts_y;
		std::vector<T> weights;
		std::vector<T> U;

				
};



#endif