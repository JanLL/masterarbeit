function [fit_data] = GN_save_fit(...
    path_root, simulation_data_struct, dsc_data_struct, index_T_dsc, ...
    optimization_data_struct, p_optim_end, optim_output)
% TODO: description!


% check if directory where we want to save our data exist
% if not, create a folder
if exist(path_root, 'dir') == 0
    disp('Creating root directory');
    mkdir(path_root);
end


fit_data = struct();


% Simulation parameters
fit_data.simulation = simulation_data_struct;

% Measurement
fit_data.measurement.dsc_data = dsc_data_struct;
fit_data.measurement.index_T_dsc = index_T_dsc;

% optimization
for fn = fieldnames(optimization_data_struct)'
   fit_data.optimization.(fn{1}) = optimization_data_struct.(fn{1});
end
fit_data.optimization.p_optim_end = p_optim_end.';
fit_data.optimization.residuum_end = optim_output.residuum_end;
fit_data.optimization.dqdp_end = optim_output.jacobian_end;

fit_data.optimization.progress_F1_norm = optim_output.progress_F1_norm;
fit_data.optimization.progress_dx_norm = optim_output.progress_dx_norm;
fit_data.optimization.progress_t_k = optim_output.progress_t_k;
fit_data.optimization.progress_NOC1 = optim_output.progress_NOC1;

fit_data.optimization.optim_duration = optim_output.optim_duration;


%%%% Save .m data file
% generic filename generation
dsc_fileSpec = dsc_data_struct.fileSpec;

datetime_cell = num2cell(int32(clock));
datetime_str = sprintf('%04i-%02i-%02i_%02i:%02i:%02i', datetime_cell{:});

heat_rate_str = strrep(num2str(dsc_data_struct.Tinfo.Tstep), '.', ',');

mass_code_str = dsc_fileSpec(11:13);

L1_str = strrep(num2str(simulation_data_struct.L1), '.', ',');
L3_str = strrep(num2str(simulation_data_struct.L3), '.', ',');


generic_fit_info_str = strcat(datetime_str, '_', mass_code_str, ...
                              '_', heat_rate_str, 'Kmin', ...
                              '_L1=', L1_str, '_L3=', L3_str);
path_fit_data_dir = strcat(path_root, generic_fit_info_str, '/');
mkdir(path_fit_data_dir);

path_data_file = strcat(path_fit_data_dir, 'fit_data.mat');
save(path_data_file, '-struct', 'fit_data');


%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%
fig = figure();
% c_p(T) plot
T_plot = 30:0.01:160;

switch optimization_data_struct.c_p_param_type
    case 'old_atan_formula'
        c_p_plot = c_p_formula(T_plot, p_optim_end(1:6));
    case 'fraser_suzuki'
        c_p_plot = c_p_fs(T_plot, p_optim_end);
    case 'gauss_linear_comb'
        c_p_plot = c_p_gauss_linear_comb(T_plot, p_optim_end);
end

plot(T_plot, c_p_plot, 'DisplayName', 'c_p Simulation');
legend('show', 'location', 'northwest');
xlabel('T [degC]');
ylabel('c_p [mJ/(mg*K]');
title(sprintf('Heat rate: %s', heat_rate_str));
path_plot_c_p = strcat(path_fit_data_dir, 'c_p(T).fig');
savefig(fig, path_plot_c_p); 

% q_pcm_in(T_ref) plot
T_ref_dsc = dsc_data_struct.data(index_T_dsc(1):index_T_dsc(2),1);
q_dsc = (dsc_data_struct.data(index_T_dsc(1):index_T_dsc(2),3) ...
      ./ dsc_data_struct.data(index_T_dsc(1):index_T_dsc(2),4)) * dsc_data_struct.mass;
  
switch simulation_data_struct.heat_rate
    case 20
        T_ref_dsc = T_ref_dsc(1:1:end);
        q_dsc = q_dsc(1:1:end);
    case 10
        T_ref_dsc = T_ref_dsc(1:2:end);
        q_dsc = q_dsc(1:2:end);
    case 5
        T_ref_dsc = T_ref_dsc(1:4:end);
        q_dsc = q_dsc(1:4:end);
    case 2.5
        T_ref_dsc = T_ref_dsc(1:8:end);
        q_dsc = q_dsc(1:8:end);
    case 1.25
        T_ref_dsc = T_ref_dsc(1:16:end);
        q_dsc = q_dsc(1:16:end);
    case 0.6
        T_ref_dsc = T_ref_dsc(1:32:end);
        q_dsc = q_dsc(1:32:end);
    case 0.3
        T_ref_dsc = T_ref_dsc(1:64:end);
        q_dsc = q_dsc(1:64:end);
    otherwise
        error('Heat rate invalid!')
end

  
q_sim = optim_output.residuum_end + q_dsc;

clf;
plot(T_ref_dsc, q_sim, 'DisplayName', 'Simulation'); hold on
plot(T_ref_dsc, q_dsc, 'DisplayName', 'Measurement');
plot(T_ref_dsc, optim_output.residuum_end, 'DisplayName', 'Residuum');
legend('show', 'location', 'northwest');
xlabel('T_{ref} [degC]');
ylabel('q_{pcm}^{in} [mW]');
title(sprintf('Heat rate: %s', heat_rate_str));
path_plot_q_pcm_in = strcat(path_fit_data_dir, 'q_pcm_in(T_ref).fig');
savefig(fig, path_plot_q_pcm_in); 

% dqdp plot
clf;
image(optim_output.jacobian_end, 'CDataMapping', 'scaled');
colorbar();
title('dq/dp')
xlabel('c_p parameters');
ylabel('T_{ref}');
set(gca,'Ydir','Normal')
path_plot_dqdp = strcat(path_fit_data_dir, 'dqdp.fig');
savefig(fig, path_plot_dqdp); 

close(fig);





end

