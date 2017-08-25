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


	private:
		int num_cntrl_pts;
		int nurbs_order;

		std::vector<double> cntrl_pts_x;
		std::vector<double> cntrl_pts_y;
		std::vector<double> weights;
		std::vector<double> U;

				
};