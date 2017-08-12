function [residuum, Jacobian] = compute_residuum_q(p_optim_free, optim_type, p_optim_estimable, p_optim_fixed, ...
                                       p_sim, q_dsc, c_p_meas, m_pcm, scalar_output, ax1, ax2)
% TODO: description!
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
c_p_const = p_sim.c_p_test_setup(1);
rho_const = p_sim.rho_test_setup(1);
lambda_const = p_sim.lambda_test_setup(1);
L1 = p_sim.L1;
N1 = p_sim.N1;
T_0 = p_sim.T_0;
T_end = p_sim.T_end;
heat_rate = p_sim.heat_rate / 60; % [K/s]



persistent T_ref T_ref_setup;

if isempty(T_ref) || isempty(T_ref_setup) || ...
   any(T_ref_setup ~= [c_p_const, rho_const, lambda_const, L1, T_0, ...
                       T_end, heat_rate])

    % Analytical
    % time grid because boundary condition T(1) = T_0 + beta * t
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

p_sim = update_c_p(p_sim, p_optim_all);

q_pcm_in = compute_q_pcm_in(p_sim, T_ref, q_dsc, m_pcm, optim_type);

residuum = q_dsc(:,2) - q_pcm_in;

% Compute Jacobian by using imaginary-Trick
% h = 0.00001;
% Jacobian = zeros(length(q_dsc(:,2)), length(p_optim_free));
% 
% for n=1:length(p_optim_free)
%     
%     p_optim_free_perturb = p_optim_free;
%     p_optim_free_perturb(n) = p_optim_free_perturb(n) + h*1i;
%     
%     p_optim_all_perturb = p_optim_all;
%     p_optim_all_perturb(p_optim_estimable) = p_optim_free_perturb;
%     
%     % update c_p evaluation functions in simulation parameter struct p_sim with
%     % new optimization parameters
%     p_sim = update_c_p(p_sim, p_optim_all_perturb);
% 
%     q_pcm_in = compute_q_pcm_in(p_sim, T_ref, q_dsc, m_pcm, optim_type);
% 
%     residuum = q_dsc(:,2) - q_pcm_in;
% 
%     Jacobian(:,n) = imag(residuum)/h;
%     
% end
% 
% figure()
% hold on;
% pcolor(Jacobian);
% shading flat % disable grid because too many lines
% colormap hsv
% colorbar;
% xlabel('Temp c_p control points [degC]');
% ylabel('Temp measurement points [degC]');


cla(ax1); % q_pcm_in plot
hold(ax1, 'on')
plot(ax1, q_dsc(:,1), q_pcm_in, 'DisplayName', 'Optimization');
plot(ax1, q_dsc(:,1), q_dsc(:,2), 'DisplayName', 'Measurement');
plot(ax1, q_dsc(:,1), residuum, 'DisplayName', 'Residuum');
legend(ax1, 'show', 'location', 'northoutside');
xlabel(ax1, 'T_{ref} [degC]');
ylabel(ax1, 'q_{pcm}^{in} [mW]');
drawnow;

cla(ax2); % c_p plot
hold(ax2, 'on')
plot(ax2, q_dsc(:,1), 1. * p_sim.eval_c_p(q_dsc(:,1)), 'DisplayName', 'Optimization'); hold on
plot(ax2, c_p_meas(:,1), c_p_meas(:,2), 'DisplayName', 'Measurement');
legend(ax2, 'show', 'location', 'northwest');
xlabel(ax2, 'T_{ref} [degC]');
ylabel(ax2, 'c_p [mJ/(mg*K]');
drawnow;

if scalar_output
    residuum = sum(residuum.^2); % compute scalar for solvers like fminsearch
end

end

