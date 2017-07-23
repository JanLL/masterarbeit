function [p_sim_mod] = update_c_p(p_sim, p_optim)
% Updates evaluation functions of c_p(T) in simulation setup struct p_sim.
%
% INPUT:
%   p_sim --> struct containing simulation setup data.
% p_optim --> optimization parameter vector containing parameters which
%             define functions c_p(T).
%
% OUTPUT:
%   p_sim_mod --> modified simulation setup data struct with updated
%                 functions to evaluate c_p(T).

if strcmp(p_sim.c_p_type, 'function_handle')
    p_sim.eval_c_p = @(T) p_sim.c_p_fcn(T, p_sim.get_param_c_p(p_optim));
    p_sim.eval_dc_p = @(T) p_sim.dc_p_fcn(T, p_sim.get_param_c_p(p_optim));
    
elseif strcmp(p_sim.c_p_type, 'B_spline')
    p_sim.c_p_bspline.knots = p_sim.get_param_c_p_knots(p_optim);
    p_sim.c_p_bspline.coefs = p_sim.get_param_c_p_coeffs(p_optim);
    %p_sim.c_p_bspline.coefs = (p_sim.get_param_c_p_coeffs(p_optim)).^2;
    % square to avoid negative coeffs when using optimizer without bounds
    p_sim.eval_c_p = @(T) spval(p_sim.c_p_bspline, T);
    
    p_sim.dc_p_bspline = fnder(p_sim.c_p_bspline);
    p_sim.eval_dc_p = @(T) spval(p_sim.dc_p_bspline, T); 
    
elseif strcmp(p_sim.c_p_type, 'NURBS')
    p_sim.c_p_nurbs.coefs(1,:) = p_sim.get_param_c_p_cntrl_pts_x(p_optim);
    p_sim.c_p_nurbs.coefs(2,:) = p_sim.get_param_c_p_cntrl_pts_y(p_optim);
    p_sim.c_p_nurbs.knots = p_sim.get_param_c_p_knots(p_optim);
    
    p_sim.dc_p_nurbs = nrbderiv(p_sim.c_p_nurbs);


    tt = 0:0.001:1;
    [c_p_curve,dc_p_curve] = ...
        nrbdeval(p_sim.c_p_nurbs, p_sim.dc_p_nurbs, tt);

    p_sim.c_p_T_table = c_p_curve(1,:);
    p_sim.c_p_table = c_p_curve(2,:);
    p_sim.dc_p_table = dc_p_curve(2,:) ./ dc_p_curve(1,:); 

    % maybe not necessary...:
    p_sim.eval_c_p = ...
        @(T) interp1(p_sim.c_p_T_table, p_sim.c_p_table, T);
    p_sim.eval_dc_p = ...
        @(T) interp1(p_sim.c_p_T_table, p_sim.dc_p_table, T);

else
    error('in update_c_p struct p_sim neither had c_p_type ''function_handle'' nor ''B_spline''.');

end

% p_sim is given as call by value, so we need to return the modified struct
p_sim_mod = deal(p_sim);

end

