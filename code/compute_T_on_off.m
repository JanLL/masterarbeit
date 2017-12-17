function [T_on, T_off] = compute_T_on_off(params, param_type)
%

if (strcmp(param_type, 'fraser_suzuki'))
    
    % from symbolic differentiation in testing2.m
    dc_p = @(T,p1,p2,p3,p4)(p1.*p2.*p3.*log(p3.*(T-p4)+1.0).*exp(-p2.*log(p3.*(T-p4)+1.0).^2).*-2.0)./(p3.*(T-p4)+1.0);

    unscaled_params = reverse_scale_params(params, param_type);
    h = unscaled_params(1);
    r = unscaled_params(2);
    wr = unscaled_params(3);
    sr = unscaled_params(4);
    z = unscaled_params(5);
    m = unscaled_params(6);
    b = unscaled_params(7);

    c1 = log(r) / log(sr)^2;
    c2 = (sr^2 - 1) / (wr*sr);

    T_infl_1 = 1/c2 * (exp((sqrt(8*c1 + 1) - 1) / (4*c1)) + c2*z - 1);
    T_infl_2 = 1/c2 * (exp(-(sqrt(8*c1 + 1) + 1) / (4*c1)) + c2*z - 1);

    m1 = dc_p(T_infl_1, h, c1, c2, z);
    m2 = dc_p(T_infl_2, h, c1, c2, z);

    b1 = c_p_fs(T_infl_1, params) - m1*T_infl_1;
    b2 = c_p_fs(T_infl_2, params) - m2*T_infl_2;

    T_on = (b - b1) / (m1 - m);
    T_off = (b - b2) / (m2 - m);
    
elseif (strcmp(param_type, 'gauss_linear_comb'))
    error('Not implemented yet!')
end




end

