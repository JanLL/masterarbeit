function [p_sim_mod] = update_c_p(p_sim, p_optim)
%

if strcmp(p_sim(1).c_p_type, 'function_handle')
    p_sim(1).eval_c_p = @(T) c_p_fcn(T, p_sim(1).get_param_c_p(p_optim));
    p_sim(1).eval_dc_p = @(T) dc_p_fcn(T, p_sim(1).get_param_c_p(p_optim));
    
elseif strcmp(p_sim(1).c_p_type, 'B_spline')
    p_sim(1).c_p_bspline.knots = p_sim(1).get_param_c_p_knots(p_optim);
    p_sim(1).c_p_bspline.coefs = p_sim(1).get_param_c_p_coeffs(p_optim);
    p_sim(1).eval_c_p = @(T) spval(p_sim(1).c_p_bspline, T);
    
    p_sim(1).dc_p_bspline = fnder(p_sim(1).c_p_bspline);
    p_sim(1).eval_dc_p = @(T) spval(p_sim(1).dc_p_bspline, T); 
    
    
    
else
    error('in update_c_p struct p_sim neither had c_p_type ''function_handle'' nor ''B_spline''.');

end

% p_sim is given as call by value, so we need to return the modified struct
p_sim_mod = deal(p_sim);

end

