



template<typename T>
void fraser_suzuki_formula(T x, T& c_p, T& dc_p, const T* params) {

	T h  = *params; params++;
	T r  = *params; params++;
	T wr = *params; params++;
	T sr = *params; params++;
	T z  = *params; params++;
	T m  = *params; params++;
	T b  = *params; params++;

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
	
	c_p_case0 = 0.01*m*x + b;
	c_p_case1 = h*exp(exp_arg) + 0.01*m*x + b;
	condassign(c_p, condition, c_p_case1, c_p_case0);

	dc_p_case0 = 0.01*m;
	dc_p_case1 = -2.*log(r)/(log_sr*log_sr) * (sr*sr - 1)/(wr*sr) * log(log_arg)/log_arg * h*exp(exp_arg) + 0.01*m;
	condassign(dc_p, condition, dc_p_case1, dc_p_case0);

	return;
}


template<typename T>
void c_p_formula(T x, T& c_p, T& dc_p, const T* params) {

	T p0 = *params; params++;
	T p1 = *params; params++;
	T p2 = *params; params++;
	T p3 = *params; params++;
	T p4 = *params; params++;
	T p5 = *params; params++;

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
void gauss_linear_comb_formula(T x, T& c_p, T& h, const T* params)
{
	const int n_gauss = 10;

	T p_ampl;
	T p_var;
	T p_offset;

	T scale_i;
	T gauss_i;

	c_p  = 0;
	h = 0;
	
	for (int i=0; i<n_gauss; ++i) {
		p_ampl   = *params; params++;
		p_var  = *params; params++;
		p_offset = *params; params++;

		scale_i = -1./(p_var);
		
		gauss_i = p_ampl * exp(scale_i * (x - p_offset)*(x - p_offset));

		c_p  += gauss_i;

	}

	T p_linear = *params; params++;
	T p_const  = *params; params++;

	c_p += 0.01*p_linear*x;
	c_p += p_const;


	return;

}
