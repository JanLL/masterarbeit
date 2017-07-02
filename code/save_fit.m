function [success] = save_fit(path_root, dsc_data, index_T_dsc, revMassNorm, ...
    p_sim, optim_solverName, optim_options, p_optim_start, p_optim_estimable, ...
    optim_con, p_optim_end, optim_output)
%

success = false;

% check if directory where we want to save our data exist
% if not, create a folder
if exist(path_root, 'dir') == 0
    disp('Creating root directory');
    mkdir(path_root);
end

% save all measurement, simulation and optimization data to comprehend what
% was done.
fit_data = struct();

fit_data.dsc_p_optim = struct();
fit_data.dsc_p_optim.dsc_data = dsc_data;
fit_data.dsc_p_optim.index_T_dsc = index_T_dsc;
fit_data.dsc_p_optim.revMassNorm = revMassNorm;

fit_data.simulation = p_sim;

fit_data.optimization = struct();
fit_data.optimization.solverName = optim_solverName;
fit_data.optimization.optim_options = optim_options;
fit_data.optimization.optim_con = optim_con;
fit_data.optimization.param_start = p_optim_start;
fit_data.optimization.estimable = p_optim_estimable;
fit_data.optimization.param_end = p_optim_end;
fit_data.optimization.output = optim_output;

path_data = strcat(path_root, 'fit_data.mat');
save(path_data, '-struct', 'fit_data');

% plot and save dU(T_ref) and c_p(T_ref) graphs

% abbreviations
heat_rate = p_sim.heat_rate / 60; % [K/min] -> [K/s]
T_0 = p_sim.T_0;
T_end = p_sim.T_end;
lambda_const = p_sim.lambda_test_setup(1);
rho_const = p_sim.rho_test_setup(1);
c_p_const = p_sim.c_p_test_setup(1);
L1 = p_sim.L1;
N1 = p_sim.N1;

% dU computation analog to compute_residuum (see there for details), 
% necessary for dU(T_ref) plot
dt = 0.05 / heat_rate; % fct evaluation every 0.05K
t = 0:dt:1/heat_rate*(T_end - T_0);
n = 100;
a = lambda_const / (c_p_const * rho_const);
T_ref = analytical_sol(L1,t,n,T_0, heat_rate, a);  

T_pcm = simulate_1d(p_sim.eval_c_p, p_sim.eval_dc_p, p_sim);

dT = T_ref(:) - T_pcm(:,N1);

k = p_sim.get_param_k(p_optim_end);
dU = polyval(k, dT);

index_T_p5 = find(T_ref(:) > p_sim.T_0 + 5, 1, 'first');

U_dsc = [dsc_data.data(index_T_dsc(1):index_T_dsc(2),1), dsc_data.data(index_T_dsc(1):index_T_dsc(2),3)];
dU_interp = interp1(T_ref(index_T_p5:end), dU(index_T_p5:end), ...
                    U_dsc(:,1), 'linear');

residuum = U_dsc(:,2) - dU_interp;

% dU(T_ref) plot
figure();
hold on;
plot(U_dsc(:,1), dU_interp, 'DisplayName', 'Optimization');
plot(U_dsc(:,1), U_dsc(:,2), 'DisplayName', 'Measurement');
plot(U_dsc(:,1), residuum, 'DisplayName', 'Residuum');
legend('show', 'location', 'northoutside');
xlabel('T_{ref} [degC]');
ylabel('\Delta U [uV]');

path_plot_dU = strcat(path_root, 'dU(T_ref).fig');
savefig(path_plot_dU);

% c_p(T_ref) plot
c_p_meas = calc_cp();

clf;
plot(U_dsc(:,1), 1. * p_sim.eval_c_p(U_dsc(:,1)), 'DisplayName', 'Optimization'); hold on
plot(c_p_meas(:,1), c_p_meas(:,2), 'DisplayName', 'Measurement');
legend('show', 'location', 'northwest');
xlabel('T_{ref} [degC]');
ylabel('c_p [mJ/(mg*K]');

path_plot_c_p = strcat(path_root, 'c_p(T_ref).fig');
savefig(path_plot_c_p);

close();

success = true;

end

