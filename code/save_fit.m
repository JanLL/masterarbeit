function [success] = save_fit(path_root, dsc_data, index_T_dsc, revMassNorm, ...
    p_sim, optim_solverName, optim_options, p_optim_start, p_optim_estimable, ...
    optim_con, p_optim_end, optim_output)
%

success = false;

% save all measurement, simulation and optimization data to comprehend what
% was done.
fit_data = struct();

fit_data.dsc_measurements = struct();
fit_data.dsc_measurements.dsc_data = dsc_data;
fit_data.dsc_measurements.dsc_data = index_T_dsc;
fit_data.dsc_measurements.dsc_data = revMassNorm;

fit_data.simulation = p_sim;

fit_data.optimization = struct();
fit_data.optimization.solverName = optim_solverName;
fit_data.optimization.optim_options = optim_options;
fit_data.optimization.optim_con = optim_con;
fit_data.optimization.param_start = p_optim_start;
fit_data.optimization.estimable = p_optim_estimable;
fit_data.optimization.param_end = p_optim_end;
fit_data.optimization.output = optim_output;

path_data = strcat(path_root, 'data.mat');
save(path_data, '-struct', 'fit_data');

% plot and save dU(T_ref) and c_p(T_ref) graphs





success = true;

end

