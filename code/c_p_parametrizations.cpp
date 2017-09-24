



template<typename T>
void fraser_suzuki_formula(T x, T& c_p, T& dc_p, T h,T r,T wr,T sr,T z, T b) {


	T c_p_case0;
	T c_p_case1;
	T dc_p_case0;
	T dc_p_case1;
	T condition;

	T log_sr = log(sr);

	// TODO: Mal mit fmin arbeiten statt komplett condition in condassign zu stecken...
	condition = z - wr*sr/(sr*sr-1) - x;

	T log_arg = (1 + (x-z)*(sr*sr-1)/(wr*sr));
	T exp_arg = -log(r)/(log_sr*log_sr) * log(log_arg)*log(log_arg);
	
	c_p_case0 = b;
	c_p_case1 = h*exp(exp_arg) + b;
	condassign(c_p, condition, c_p_case1, c_p_case0);

	dc_p_case0 = 0.;
	dc_p_case1 = -2.*log(r)/(log_sr*log_sr) * (sr*sr - 1)/(wr*sr) * log(log_arg)/log_arg * h*exp(exp_arg);
	condassign(dc_p, condition, dc_p_case1, dc_p_case0);

	return;
}


template<typename T>
void c_p_formula(T x, T& c_p, T& dc_p, T p0, T p1, T p2, T p3, T p4, T p5) {

	T atan_arg = -p3*(x-p0);
	T atan_value = atan(atan_arg);
	T exp_arg = -p2*(x-p0)*(x-p0);
	T exp_value = exp(exp_arg);

	c_p = (atan_value + M_PI/2.)*(p1*exp_value) + p4*x + p5;

	dc_p = -(p1*p3*exp_value)/(1+atan_value*atan_value) 
	       -2*p1*p2*(x-p0)*exp_value * (atan_value + M_PI/2.) + p4;

	return;
}



template<typename T>
void gauss_linear_comb_formula(T x, T& c_p, T& dc_p, 
		T ampl_0, T sigma_0, T shift_0,
		T ampl_1, T sigma_1, T shift_1,
		T ampl_2, T sigma_2, T shift_2,
		T ampl_3, T sigma_3, T shift_3,
		T ampl_4, T sigma_4, T shift_4,
		T ampl_5, T sigma_5, T shift_5,
		T ampl_6, T sigma_6, T shift_6,
		T ampl_7, T sigma_7, T shift_7,
		T ampl_8, T sigma_8, T shift_8,
		T ampl_9, T sigma_9, T shift_9,
		T m, T b)
{

	T scale_0 = -1./(sigma_0*sigma_0);
	T scale_1 = -1./(sigma_1*sigma_1);
	T scale_2 = -1./(sigma_2*sigma_2);
	T scale_3 = -1./(sigma_3*sigma_3);
	T scale_4 = -1./(sigma_4*sigma_4); 
	T scale_5 = -1./(sigma_5*sigma_5); 
	T scale_6 = -1./(sigma_6*sigma_6); 
	T scale_7 = -1./(sigma_7*sigma_7); 
	T scale_8 = -1./(sigma_8*sigma_8); 
	T scale_9 = -1./(sigma_9*sigma_9); 

	T c_p_0 = ampl_0 * exp(scale_0 * (x - shift_0)*(x - shift_0));
	T c_p_1 = ampl_1 * exp(scale_1 * (x - shift_1)*(x - shift_1));
	T c_p_2 = ampl_2 * exp(scale_2 * (x - shift_2)*(x - shift_2));
	T c_p_3 = ampl_3 * exp(scale_3 * (x - shift_3)*(x - shift_3));
	T c_p_4 = ampl_4 * exp(scale_4 * (x - shift_4)*(x - shift_4));
	T c_p_5 = ampl_5 * exp(scale_5 * (x - shift_5)*(x - shift_5));
	T c_p_6 = ampl_6 * exp(scale_6 * (x - shift_6)*(x - shift_6));
	T c_p_7 = ampl_7 * exp(scale_7 * (x - shift_7)*(x - shift_7));
	T c_p_8 = ampl_8 * exp(scale_8 * (x - shift_8)*(x - shift_8));
	T c_p_9 = ampl_9 * exp(scale_9 * (x - shift_9)*(x - shift_9));

	c_p = c_p_0 + c_p_1 + c_p_2 + c_p_3 + c_p_4 + c_p_5 + c_p_6 + c_p_7 + c_p_8 + c_p_9 + 0.01*m*x + b;

	dc_p = 	 2*scale_0 * (x-shift_0) * c_p_0
		   + 2*scale_1 * (x-shift_1) * c_p_1
		   + 2*scale_2 * (x-shift_2) * c_p_2
		   + 2*scale_3 * (x-shift_3) * c_p_3
		   + 2*scale_4 * (x-shift_4) * c_p_4
		   + 2*scale_5 * (x-shift_5) * c_p_5
		   + 2*scale_6 * (x-shift_6) * c_p_6
		   + 2*scale_7 * (x-shift_7) * c_p_7
		   + 2*scale_8 * (x-shift_8) * c_p_8
		   + 2*scale_9 * (x-shift_9) * c_p_9
		   + 0.01*m;

	return;

}
