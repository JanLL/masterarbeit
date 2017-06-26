function [residuum] = compute_residuum(p_optim_free, p_optim_estimable, p_optim_fixed, ...
                                       p_sim, U_dsc, c_p_meas, ax1, ax2)
% Given DSC measurements U_dsc(T_ref) and simulation setup. From the
% simulation we also gain a voltage U_sim(T_ref). From these two quantities
% we can compute a residuum between measurements and simulation.
%
% INPUT:
%  p_optim_free --> free optimization variables:
%                    1) parameters which define function c_p(T)
%                    2) parameters k which define function U(T)
% p_optim_estimable --> boolean array to specify which optimization
%                       variables are free (true) or fixed(false)
% p_optim_fixed --> values of fixed optimization parameters. Together with 
%                   p_optim_estimable and p_optim one can reconstruct the
%                   whole (free and fixed) optimization variable vector.
%      p_sim --> struct with simulation parameters, available via the
%                function get_default_sim_params()
%      U_dsc --> array where U_dsc(:,1) contains temperatures T_ref and 
%                U_dsc(:,2) corresponding voltages from DSC measurement.
%   c_p_meas --> array where c_p_meas(:,1) contains temperatures T_ref and
%                c_p_meas(:,2) corresponding calculated c_p values from DSC
%                measurement.
%        ax1 --> axis of dU(T_ref) plot.
%        ax2 --> axis of c_p(T_ref) plot.
%
% OUTPUT:
%   residuum --> residuum between voltage from measurement and simulation.


% abbreviations
c_p_const = p_sim(1).c_p_test_setup(1);
rho_const = p_sim(1).rho_test_setup(1);
lambda_const = p_sim(1).lambda_test_setup(1);
L1 = p_sim(1).L1;
N1 = p_sim(1).N1;
T_0 = p_sim(1).T_0;
T_end = p_sim(1).T_end;
heat_rate = p_sim(1).heat_rate / 60; % [K/s]



persistent T_ref T_ref_setup;

if isempty(T_ref) || isempty(T_ref_setup) || ...
   any(T_ref_setup ~= [c_p_const, rho_const, lambda_const, L1, T_0, ...
                       T_end, heat_rate])

    % Analytical
    dt = 0.05 / heat_rate; % fct evaluation every 0.05K
    t = 0:dt:1/heat_rate*(T_end - T_0);
    n = 100;
    a = lambda_const / (c_p_const * rho_const);
    
    T_ref = analytical_sol(L1,t,n,T_0, heat_rate, a);  
    T_ref_setup = [c_p_const, rho_const, lambda_const, L1, T_0, ...
                   T_end, heat_rate];
end

% build vector of all (free and fixed) optimization parameters
p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim_free;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

% update c_p evaluation functions in simulation parameter struct p_sim with
% new optimization parameters
p_sim = update_c_p(p_sim, p_optim_all);

% compute temperature difference between pcm and reference
%T_pcm = simulate_1d(eval_c_p_expl, eval_dc_p_expl, p_sim(1));
T_pcm = simulate_1d(p_sim(1).eval_c_p, p_sim(1).eval_dc_p, p_sim(1));


dT = T_ref(:) - T_pcm(:,N1);

% convert temperatue to voltage difference with linear factor k
k = p_sim(1).get_param_k(p_optim_all);

dU = k(1) + k(2) * dT + k(3) * dT.^2;

% In the first few seconds T_ref(:,N1) stays constant till the heat of the
% oven reaches it. Interpolation can't handle non-unique domain of
% definition.
index_T_p5 = find(T_ref(:) > p_sim(1).T_0 + 5, 1, 'first');
dU_interp = interp1(T_ref(index_T_p5:end), dU(index_T_p5:end), ...
                    U_dsc(:,1), 'linear');

% Note: Interpolation necessary because T_ref doesnt increase linearly, at
% least at the beginning, i.e. we dont start at T_0, there's a delay of the
% Constantan.

residuum = U_dsc(:,2) - dU_interp;

cla(ax1);
hold(ax1, 'on')
plot(ax1, U_dsc(:,1), dU_interp, 'DisplayName', 'Optimization');
plot(ax1, U_dsc(:,1), U_dsc(:,2), 'DisplayName', 'Measurement');
plot(ax1, U_dsc(:,1), residuum, 'DisplayName', 'Residuum');
legend(ax1, 'show', 'location', 'northoutside');
xlabel(ax1, 'T_{ref}');
ylabel(ax1, '\Delta U');
drawnow;

cla(ax2);
hold(ax2, 'on')
plot(ax2, U_dsc(:,1), p_sim(1).eval_c_p(U_dsc(:,1)), 'DisplayName', 'Optimization'); hold on
plot(ax2, c_p_meas(:,1), c_p_meas(:,2), 'DisplayName', 'Measurement');
legend(ax2, 'show', 'location', 'northwest');
xlabel(ax2, 'T_{ref}');
ylabel(ax2, 'c_p');
drawnow;

% % print('/home/argo/masterarbeit/simulationen-data/delta_U_optimized_lambda-wikiX8','-dpng');

end

