function [fit_data] = save_fit(path_root, dsc_data_struct, index_T_dsc, revMassNorm, ...
    p_sim, optim_solverName, optim_options, p_optim_start, p_optim_estimable, ...
    optim_con, p_optim_end, optim_output)
%


% check if directory where we want to save our data exist
% if not, create a folder
if exist(path_root, 'dir') == 0
    disp('Creating root directory');
    mkdir(path_root);
end

% save all measurement, simulation and optimization data to comprehend what
% was done.
fit_data = struct();

fit_data.measurements = struct();
fit_data.measurements.dsc_data_struct = dsc_data_struct;
fit_data.measurements.index_T_dsc = index_T_dsc;
fit_data.measurements.revMassNorm = revMassNorm;

fit_data.simulation = p_sim;

fit_data.optimization = struct();
fit_data.optimization.solverName = optim_solverName;
fit_data.optimization.optim_options = optim_options;
fit_data.optimization.optim_con = optim_con;
fit_data.optimization.param_start = p_optim_start;
fit_data.optimization.estimable = p_optim_estimable;
fit_data.optimization.param_end = p_optim_end;
fit_data.optimization.output = optim_output;

% generic filename generation
dsc_fileSpec = dsc_data_struct.fileSpec;

datetime_cell = num2cell(int32(clock));
datetime_str = sprintf('%04i-%02i-%02i_%02i:%02i:%02i', datetime_cell{:});

heat_rate_str = num2str(p_sim.heat_rate);
heat_rate_str = strrep(heat_rate_str, '.', ',');

mass_code_str = dsc_fileSpec(11:13);

generic_fit_info_str = strcat(datetime_str, '_', mass_code_str, '_', heat_rate_str, ...
                   'Kmin_', optim_solverName);
path_fit_data_dir = strcat(path_root, generic_fit_info_str, '/');
mkdir(path_fit_data_dir);

path_data_file = strcat(path_fit_data_dir, 'fit_data.mat');
save(path_data_file, '-struct', 'fit_data');


% plot and save dU(T_ref) and c_p(T_ref) graphs
dU = compute_dU(fit_data);
residuum = dU(:,2) - dU(:,3);

% dU(T_ref) plot
figure();
hold on;
plot(dU(:,1), dU(:,3), 'DisplayName', 'Optimization');
plot(dU(:,1), dU(:,2), 'DisplayName', 'Measurement');
plot(dU(:,1), residuum, 'DisplayName', 'Residuum');
legend('show', 'location', 'northoutside');
xlabel('T_{ref} [degC]');
ylabel('\Delta U [uV]');

path_plot_dU = strcat(path_fit_data_dir, 'dU(T_ref).fig');
savefig(path_plot_dU);

% c_p(T_ref) plot

c_p_meas = calc_cp(dsc_data_struct);
c_p_optim = p_sim.eval_c_p(dU(:,1));

factor = max(c_p_meas(:,2)) / max(c_p_optim);

clf;
plot(dU(:,1), c_p_optim, 'DisplayName', 'Optimization'); hold on
plot(c_p_meas(:,1), c_p_meas(:,2), 'DisplayName', 'Measurement');
plot(dU(:,1), factor * c_p_optim, 'DisplayName', sprintf('Optimization X %g', factor)); hold on
legend('show', 'location', 'northwest');
xlabel('T_{ref} [degC]');
ylabel('c_p [mJ/(mg*K]');

path_plot_c_p = strcat(path_fit_data_dir, 'c_p(T_ref).fig');
savefig(path_plot_c_p);

close();


end

