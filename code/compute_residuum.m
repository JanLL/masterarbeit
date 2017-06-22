function [residuum] = compute_residuum(p_optim, p_optim_estimable, p_optim_fixed, ...
                                       p_sim, U_dsc, T_meas, c_p_meas, T_ref_meas, ax1, ax2)
% TODO: function description!

%p_optim

% abbreviations
c_p_const = p_sim(1).c_p_test_setup(1);
rho_const = p_sim(1).rho_test_setup(1);
lambda_const = p_sim(1).lambda_test_setup(1);
L1 = p_sim(1).L1;
N1 = p_sim(1).N1;
T_0 = p_sim(1).T_0;
T_end = p_sim(1).T_end;
heat_rate = p_sim(1).heat_rate;



persistent T_ref T_ref_setup;

if isempty(T_ref) || isempty(T_ref_setup) || ...
   any(T_ref_setup ~= [c_p_const, rho_const, lambda_const, L1, N1, T_0, ...
                       T_end, heat_rate])
       
    T_ref = simulate_1d(p_sim(1).eval_c_p, p_sim(1).eval_dc_p, p_sim(2));
    T_ref_setup = [c_p_const, rho_const, lambda_const, L1, N1, T_0, ...
                   T_end, heat_rate];
end

% build vector of all (free and fixed) optimization parameters
p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

% update c_p evaluation functions in simulation parameter struct p_sim with
% new optimization parameters
p_sim = update_c_p(p_sim, p_optim_all);

% compute temperature difference between pcm and reference
%T_pcm = simulate_1d(eval_c_p_expl, eval_dc_p_expl, p_sim(1));
T_pcm = simulate_1d(p_sim(1).eval_c_p, p_sim(1).eval_dc_p, p_sim(1));

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

cla(ax1);
hold(ax1, 'on')
plot(ax1, T_ref_meas, dU_interp, 'DisplayName', 'Optimization');
plot(ax1, T_ref_meas, U_dsc, 'DisplayName', 'Measurement');
plot(ax1, T_ref_meas, residuum, 'DisplayName', 'Residuum');
legend(ax1, 'show', 'location', 'northwest');
xlabel(ax1, 'T_{ref}');
ylabel(ax1, '\Delta U');
drawnow;

cla(ax2);
hold(ax2, 'on')
plot(ax2, T_ref_meas, p_sim(1).eval_c_p(T_ref_meas), 'DisplayName', 'Optimization'); hold on
plot(ax2, T_meas, c_p_meas, 'DisplayName', 'Measurement');
legend(ax2, 'show', 'location', 'northwest');
xlabel(ax2, 'T_{ref}');
ylabel(ax2, 'c_p');
drawnow;

% % print('/home/argo/masterarbeit/simulationen-data/delta_U_optimized_lambda-wikiX8','-dpng');

end

