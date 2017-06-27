tic;

T0 = 30.;
heat_rate = 10. / 60; % K/s
L = 25.;

c_p = 0.41;
rho = 8.9;
lambda = 23.;

a = lambda / (c_p * rho);

x = [0, L];
n = 100;

T_ref_max = 200';

% Solve non-linear system to get time t where T_ref (temp. at crucible)
% reaches T_ref_max.
% F = @(t) analytical_sol(x,t,n,T0, heat_rate, a) - T_ref_max;
% 
% fsolve_options = optimoptions('fsolve','Display','none');
% t_max = fsolve(F, t, fsolve_options);
% 
% dt = 0.05 / heat_rate; % fct evaluation every 0.05K
% t = 0:dt:t_max;

dt = 0.05 / heat_rate; % fct evaluation every 0.05K
t = 0:dt:1/heat_rate*(T_ref_max - T0);

T_ref_ana = analytical_sol(x,t,n,T0, heat_rate, a);

dT = T_ref_ana(:,1) - T_ref_ana(:,2);
toc

tic;
% Numerical Solution
% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N1 = 1250;
N2 = 0;
N3 = 0;

T_0 = 30;
T_end = 200;

heat_rate = 10.; % K/min

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N1', N1, 'N2', N2, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate};
p_sim = get_param_sim(common_args{:});


c_p_params = [144.0009 - 15., ...
                    4.1036 * 5., ...
                    0.0039 + 0.1, ...
                    1.4217 * 0., ...
                    0.0078, ...
                    1.5325];

p_sim = update_c_p(p_sim, c_p_params);

                
T_ref_num = simulate_1d(p_sim.eval_c_p, p_sim.eval_dc_p, p_sim);
toc

% T_rel_err = T_ref_ana ./ T_ref_num(:,end);
% T_abs_err = abs(T_ref_ana - T_ref_num(:,end));
% %plot(T_ref_ana, T_rel_err, '.'); hold on
% plot(T_ref_ana, T_abs_err, '.', 'DisplayName', sprintf('N1=%i', N3*L1)); hold on
% xlabel('T_{ref}(analytical)')
% ylabel('|T_{ref}(analytical) - T_{ref}(numerical)|')
% title('Absolute error with n_{ana}=100');

%legend('show', 'location', 'northoutside')


plot(T_ref_ana(:,2), dT, '--'); hold on
plot(T_ref_num(:,end), T_ref_num(:,1) - T_ref_num(:,end), '--');



