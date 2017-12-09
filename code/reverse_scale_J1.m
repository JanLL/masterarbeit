function [unscaled_J1] = reverse_scale_J1(scaled_J1, param_type, p_optim_estimable)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n_mp = size(scaled_J1,1);

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
    scaling_gausse_matrix = repmat(scaling_gausse(p_optim_estimable),n_mp,1);
    scaling_gausse_matrix = 1 ./ scaling_gausse_matrix;

    unscaled_J1 = scaled_J1 .* scaling_gausse_matrix;
    
elseif (strcmp(param_type, 'fraser_suzuki'))
    error('todo!!')
    
end





end

