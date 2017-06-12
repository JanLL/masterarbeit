function [residuum] = compute_residuum(p_optim, p_optim_estimable, p_optim_fixed, ...
                                       p_sim, U_dsc, T_ref_meas)
% TODO: function description!


% abbreviations
N1 = p_sim(1).N1;

% build vector of all (free and fixed) optimization parameters
p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

% build new function handles of c_p and derivative with explicit parameter
% values from p_optim.
c_p_params = p_sim(1).get_param_c_p(p_optim_all);
eval_c_p_expl = @(T) p_sim(1).eval_c_p(T, c_p_params);
eval_dc_p_expl = @(T) p_sim(1).eval_dc_p(T, c_p_params);

% compute temperature difference between pcm and reference
T_pcm = simulate_1d(eval_c_p_expl, eval_dc_p_expl, p_sim(1));
T_ref = simulate_1d(eval_c_p_expl, eval_dc_p_expl, p_sim(2));

dT = T_ref(:,N1) - T_pcm(:,N1);

% convert temperatue to voltage difference with linear factor k
k = p_sim(1).get_param_k(p_optim_all);

dU = k * dT;

% In the first few seconds T_ref(:,N1) stays constant till the heat of the
% oven reaches it. Interpolation can't handle non-unique domain of
% definition.
index_T_p5 = find(T_ref(:,N1) > p_sim(1).T_0 + 5, 1, 'first');
dU_interp = interp1(T_ref(index_T_p5:end,N1), dU(index_T_p5:end), ...
                    T_ref_meas, 'linear');

% Note: Interpolation necessary because T_ref doesnt increase linearly, at
% least at the beginning, i.e. we dont start at T_0, there's a delay of the
% Constantan.

residuum = U_dsc - dU_interp;

% figure()
% plot(T_ref_meas, dU_interp, 'DisplayName', 'Simulation'); hold on
% plot(T_ref_meas, U_dsc, 'DisplayName', 'Measurements'); hold on
% plot(T_ref_meas, residuum, 'DisplayName', 'Residuum'); hold on
% legend('show', 'location', 'northwest');
% xlabel('T_{ref}');
% ylabel('\Delta U');

% % print('/home/argo/masterarbeit/simulationen-data/delta_U_optimized_lambda-wikiX8','-dpng');

end

