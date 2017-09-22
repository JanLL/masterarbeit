function [fit_data] = save_fit_mex(...
    path_root, simulation_data_struct, dsc_data_struct, index_T_dsc, ...
    optimization_data_struct, p_optim_end, optim_output, optim_duration, ...
    residuum_end, Jacobian_end)
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
fit_data.optimization.p_optim_end = p_optim_end;
fit_data.optimization.optim_output = optim_output;
fit_data.optimization.optim_duration = optim_duration;
fit_data.optimization.dqdp_end = Jacobian_end;

%%%% Save .m data file
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

%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%
figure();
% c_p(T)
T_plot = 30:0.01:160;
%c_p_plot = c_p_fs(T_plot, p_optim_all);
%c_p_plot = c_p_formula(T_plot, p_optim_all(1:6));
c_p_plot = c_p_gauss_linear_comb(T_plot, p_optim_all);
plot(aT_plot, c_p_plot, 'DisplayName', 'c_p Simulation');
legend('show', 'location', 'northoutside');
xlabel('T [degC]');
ylabel('c_p [mJ/(mg*K]');
path_plot_c_p = strcat(path_fit_data_dir, 'c_p(T).fig');
savefig(fig, path_plot_c_p); 

% q_pcm_in(T_ref)
dsc = simulation_data_struct;
T_ref_dsc = dsc.data(index_T_dsc(1):index_T_dsc(2),1);
q_dsc = (dsc.data(index_T_dsc(1):index_T_dsc(2),3) ...
      ./ dsc.data(index_T_dsc(1):index_T_dsc(2),4)) * m_pcm;
q_sim = residuum_end + q_dsc;

clf;
plot(T_ref_dsc, q_sim, 'DisplayName', 'Simulation');
plot(T_ref_dsc, q_dsc, 'DisplayName', 'Measurement');
plot(T_ref_dsc, residuum_end, 'DisplayName', 'Residuum');
legend('show', 'location', 'northoutside');
xlabel('T_{ref} [degC]');
ylabel('q_{pcm}^{in} [mW]');
path_plot_q_pcm_in = strcat(path_fit_data_dir, 'q_pcm_in(T_ref).fig');
savefig(fig, path_plot_q_pcm_in); 

% dqdp 
clf;
image(Jacobian_end, 'CDataMapping', 'scaled');
colorbar();
title('dq/dp')
xlabel('c_p parameters');
ylabel('T_{ref}');
path_plot_dqdp = strcat(path_fit_data_dir, 'dqdp.fig');
savefig(fig, path_plot_dqdp); 

close();



end

