function [residuum] = compute_residuum(p_optim, p_sim)
%



% compute derivative (symbolically) of c_p w.r.t. temperature T
persistent dc_p;

if isempty(dc_p)
    syms T;
    p = sym('p', [6 1]);
    dc_p = matlabFunction(diff(c_p_formula(T, p), T), 'Vars', [T;p]);
end

% build new function handles of c_p and derivative with explicit parameter
% values from p_optim.
c_p_params = get_c_p_params(p_optim);
eval_c_p = @(T)c_p_formula(T, c_p_params);
c_p_params_cell = num2cell(c_p_params);
eval_dc_p = @(T) dc_p(T,c_p_params_cell{1:end});





end

