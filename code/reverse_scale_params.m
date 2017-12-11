function [unscaled_params] = reverse_scale_params(scaled_params, param_type)
%

% Scaling values Gaussians
scaling_ampl   = [2.3314, 13.3911, -2.6634, 1, 1, 1, 1, 1, 1, 1];
scaling_var    = [473.9, 24.51, 121.8, 30, 30, 30, 30, 30, 30, 30];
scaling_offset = [134.38, 127.9, 145.7, 130, 130, 130, 130, 130, 130, 130];
scaling_linear = 0.01;
scaling_const  = 1.;

scaling_gausse_tmp = [scaling_ampl; scaling_var; scaling_offset];
scaling_gausse = [reshape(scaling_gausse_tmp,1,[]), scaling_linear, scaling_const];

% Scaling values Fraser-Suzuki
scale_h  = 14.;
scale_r  = 2.;
scale_wr = 10.7;
scale_sr = 0.705;
scale_z  = 129.;
scale_m  = 0.00789;
scale_b  = 1.69;

scaling_fs = [scale_h, scale_r, scale_wr, scale_sr, scale_z, scale_m, scale_b];

if (strcmp(param_type, 'gauss_linear_comb')) 
    unscaled_params = scaling_gausse.' .* scaled_params;
elseif (strcmp(param_type, 'fraser_suzuki'))
    unscaled_params = scaling_fs.' .* scaled_params;
end




end

