function [fit_data] = save_fit(path_root, dsc_data_struct, index_T_dsc, revMassNorm, ...
    p_sim, optim_type, optim_solverName, optim_options, p_optim_start, p_optim_estimable, ...
    optim_con, p_optim_end, optim_output, optim_jac_output)
% Saves the main fit results in form of the graphs c_p(T) and dU(T_ref)
% with the optimized parameters.
% Additionally saves all necessary data to be able to reproduce this fit.
%
% INPUT: 
%       path_root --> Directory where all fits are saved. A single fit is
%                     then saved in a special directory with generic name.
% dsc_data_struct --> struct with dsc measurement information.
%     index_T_dsc --> min/max index pair of used domain of T_ref for
%                     residuum computation.
%     revMassNorm --> logical value whether mass normalization U_dsc is 
%                     reversed.
%           p_sim --> simulation parameter struct.
% optim_solverName -> name of optimization solver used.
%   optim_options --> options of optimization solver.
%   p_optim_start --> start optimization parameter values.
% p_optim_estimable -> logical array which optimization parameters are
%                      fixed (false) and free to optimize (true).
%       optim_con --> cell array with optimization constraints.
%     p_optim_end --> array with optimized optimization variables.
%    optim_output --> information output struct TODO: FNISH!!


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
fit_data.optimization.optim_type = optim_type;
fit_data.optimization.solverName = optim_solverName;
fit_data.optimization.optim_options = optim_options;
fit_data.optimization.optim_con = optim_con;
fit_data.optimization.param_start = p_optim_start;
fit_data.optimization.estimable = p_optim_estimable;
fit_data.optimization.param_end = p_optim_end;
fit_data.optimization.output = optim_output;
fit_data.optimization.jac_output = optim_jac_output;


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

if strcmp(optim_type(1:end), 'dU')
    % plot and save dU(T_ref)
    fit_quantity = compute_dU(fit_data);
    residuum = fit_quantity(:,2) - fit_quantity(:,3);

    figure();
    hold on;
    plot(fit_quantity(:,1), fit_quantity(:,3), 'DisplayName', 'Optimization');
    plot(fit_quantity(:,1), fit_quantity(:,2), 'DisplayName', 'Measurement');
    plot(fit_quantity(:,1), residuum, 'DisplayName', 'Residuum');
    legend('show', 'location', 'northoutside');
    xlabel('T_{ref} [degC]');
    ylabel('\Delta U [uV]');

    path_plot_dU = strcat(path_fit_data_dir, 'dU(T_ref).fig');
    savefig(path_plot_dU);
    
elseif strcmp(optim_type(1:end-2), 'heat_flux')
    fit_quantity = compute_q_pcm_in(fit_data);
    residuum = fit_quantity(:,2) - fit_quantity(:,3);

    figure()
    hold on;
    plot(fit_quantity(:,1), fit_quantity(:,3), 'DisplayName', 'Optimization');
    plot(fit_quantity(:,1), fit_quantity(:,2), 'DisplayName', 'Measurement');
    plot(fit_quantity(:,1), residuum, 'DisplayName', 'Residuum');
    legend('show', 'location', 'northoutside');
    xlabel('T_{ref} [degC]');
    ylabel('q_{pcm}^{in} [mW]');
    
    path_plot_q_pcm_in = strcat(path_fit_data_dir, 'q_pcm_in(T_ref).fig');
    savefig(path_plot_q_pcm_in);    
end

% c_p(T_ref) plot
c_p_meas = calc_cp(dsc_data_struct);
c_p_optim = p_sim.eval_c_p(fit_quantity(:,1));

%factor = max(c_p_meas(:,2)) / max(c_p_optim); % old
index_c_p_meas_1 = find(c_p_meas(:,1) > 30, 1, 'first');
c_p_meas_1 = c_p_meas(index_c_p_meas_1,2);
index_c_p_meas_2 = find(c_p_meas(index_c_p_meas_1:end,2) > c_p_meas_1*1.15, 1, 'first') + index_c_p_meas_1;
mean_c_p_meas = mean(c_p_meas(index_c_p_meas_1:index_c_p_meas_2,2));

index_c_p_optim_1 = find(fit_quantity(:,1) > 30, 1, 'first');
index_c_p_optim_2 = find(fit_quantity(:,1) > c_p_meas(index_c_p_meas_2,1),1,'first');
mean_c_p_optim = mean(c_p_optim(index_c_p_optim_1:index_c_p_optim_2));

factor = mean_c_p_meas / mean_c_p_optim;

clf;
plot(fit_quantity(:,1), c_p_optim, 'DisplayName', 'Optimization'); hold on
plot(c_p_meas(:,1), c_p_meas(:,2), 'DisplayName', 'Measurement');
%plot(fit_quantity(:,1), factor * c_p_optim, 'DisplayName', sprintf('Optimization X %g', factor)); hold on
legend('show', 'location', 'northwest');
xlabel('T_{ref} [degC]');
ylabel('c_p [mJ/(mg*K]');

path_plot_c_p = strcat(path_fit_data_dir, 'c_p(T_ref).fig');
savefig(path_plot_c_p);

close();

% save jacobian if lsqnonlin used (currently just works for NURBS)
if strcmp(optim_solverName, 'lsqnonlin') && strcmp(p_sim.c_p_type, 'NURBS')
    U_dsc = [dsc_data_struct.data(index_T_dsc(1):index_T_dsc(2),1), ...
             dsc_data_struct.data(index_T_dsc(1):index_T_dsc(2),3)];
    num_cntrl_pts = p_sim.c_p_params_num(1);
    cntrl_pts_x = p_optim_end(1:num_cntrl_pts);
    
    figure()
    hold on;
    pcolor(cntrl_pts_x, U_dsc(:,1), optim_jac_output);
    shading flat % disable grid because too many lines
    colormap hsv
    colorbar;
    xlabel('Temp c_p control points [degC]');
    ylabel('Temp measurement points [degC]');
    
    path_plot_jac_output = strcat(path_fit_data_dir, 'jac_output.fig');
    savefig(path_plot_jac_output);
    close();
end


end

