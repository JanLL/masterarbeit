



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
