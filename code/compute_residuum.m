function [residuum] = compute_residuum(p_optim, p_sim, U_dsc, T_ref_meas)
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
c_p_params = get_param_c_p(p_optim);
eval_c_p = @(T)c_p_formula(T, c_p_params);
c_p_params_cell = num2cell(c_p_params);
eval_dc_p = @(T) dc_p(T,c_p_params_cell{1:end});

% compute temperature difference between pcm and reference

T_pcm = simulate_1d(eval_c_p, eval_dc_p, p_sim(1));
T_ref = simulate_1d(eval_c_p, eval_dc_p, p_sim(2));

N1 = p_sim(1).N1;
dT = T_ref(:,N1) - T_pcm(:,N1);

% convert temperatue to voltage difference with linear factor k
k = get_param_k(p_optim);

dU = k * dT;

tic;
dU_interp = interp1(T_ref(:,N1), dU, T_ref_meas, 'linear');
toc

% Note: Interpolation necessary because T_ref doesnt increase linearly, at
% least at the beginning, i.e. we dont start at T_0, theres a delay of the
% Constantan.

residuum = U_dsc - dU_interp;

plot(T_ref_meas, dU_interp); hold on
plot(T_ref_meas, U_dsc); hold on
plot(T_ref_meas, residuum); hold on


end

