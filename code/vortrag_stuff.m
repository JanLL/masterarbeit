%% Fuer alle Heizraten c_p(T) und heat flux seperat plotten als .png
%  c_p erstellen fuer Latex


% fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-08_22:22:31_407_L1=40_L3=0,1_N1=300_N3=50_5Gaussians/';
% fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-09_18:33:20_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS/';

fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-19_20:27:59_407_L1=40_L3=0.1_N1=300_N3=50_5Gaussians_used/';
% fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-20_14:25:10_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS_used/';


file_list = dir(fit_dir);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];


fig1 = figure(1); clf; ax1 = gca; hold on;
set(fig1, 'Units', 'normalized', 'OuterPosition', [0., 0., 0.8, 1.]); 
fig2 = figure(2); clf; ax2 = gca; hold on;
set(fig2, 'Units', 'normalized', 'OuterPosition', [0., 0., 0.8, 1.]); 

for j=1:length(nameSubDirs)
    
    cla(ax1);
    cla(ax2);
    
    fit_data = load([fit_dir, nameSubDirs{j}, '/fit_data.mat']);
    p_optim_all = fit_data.optimization.p_optim_end;
    

    
    % c_p(T) plot
    T_plot = 30:0.01:160;
    switch fit_data.optimization.c_p_param_type
        case 'old_atan_formula'
            c_p_plot = c_p_formula(T_plot, p_optim_all(1:6));
        case 'fraser_suzuki'
            c_p_plot = c_p_fs(T_plot, p_optim_all);
        case 'gauss_linear_comb'
            c_p_plot = c_p_gauss_linear_comb(T_plot, p_optim_all);
    end    
    
    plot(ax1, T_plot, c_p_plot, 'DisplayName', 'c_p(T)', 'Linewidth', 2.)
    legend(ax1, 'show', 'location', 'northwest');
    xlabel(ax1, 'T [°C]');
    ylabel(ax1, 'c_p [mJ/(mg*K)]');
    set(ax1,'FontSize',30);
    xlim(ax1, [30, 160]);
    title(ax1, '')

    % heat flux plot
    index_T_dsc = fit_data.measurement.index_T_dsc;
    meas_data = fit_data.measurement.dsc_data.data;
    T_ref_dsc = meas_data(index_T_dsc(1):index_T_dsc(2), 1);
    
    q_res = fit_data.optimization.residuum_end;
    q_meas = fit_data.measurement.dsc_data.mass * ...
             meas_data(index_T_dsc(1):index_T_dsc(2),3) ./ ...
             meas_data(index_T_dsc(1):index_T_dsc(2),4);
         
%     switch fit_data.simulation.heat_rate
%         case 20
%             T_ref_dsc = T_ref_dsc(1:1:end);
%             q_meas = q_meas(1:1:end);
%         case 10
%             T_ref_dsc = T_ref_dsc(1:2:end);
%             q_meas = q_meas(1:2:end);
%         case 5
%             T_ref_dsc = T_ref_dsc(1:4:end);
%             q_meas = q_meas(1:4:end);
%         case 2.5
%             T_ref_dsc = T_ref_dsc(1:8:end);
%             q_meas = q_meas(1:8:end);
%         case 1.25
%             T_ref_dsc = T_ref_dsc(1:16:end);
%             q_meas = q_meas(1:16:end);
%         case 0.6
%             T_ref_dsc = T_ref_dsc(1:32:end);
%             q_meas = q_meas(1:32:end);
%         case 0.3
%             T_ref_dsc = T_ref_dsc(1:64:end);
%             q_meas = q_meas(1:64:end);
%         otherwise
%             error('Heat rate invalid!')
%     end
         
         
    q_sim = q_res + q_meas;
    
    plot(ax2, T_ref_dsc, q_sim, 'DisplayName', 'Simulation', 'Linewidth', 2.); hold on
    plot(ax2, T_ref_dsc, q_meas, 'DisplayName', 'Measurement', ...
         'Linestyle', '--', 'Linewidth', 2.);
    plot(ax2, T_ref_dsc, q_res, 'DisplayName', 'Residuum', 'Linewidth', 2.);
    legend(ax2, 'show', 'location', 'northwest');
    xlabel(ax2, 'T_{ref} [°C]');
    ylabel(ax2, '\Phi_q^{pcm,in} [mW]');
    set(ax2,'FontSize',30);
    xlim(ax2, [T_ref_dsc(1), T_ref_dsc(end)]);
    title(ax2, '')
    
    print(fig1, [fit_dir, nameSubDirs{j}, '/c_p_vortrag'], '-dpng', '-r200');
    print(fig2, [fit_dir, nameSubDirs{j}, '/heat_flux_vortrag'], '-dpng', '-r200');
    
    
end

close(1);
close(2);



%% Plot heat flux measurements and c_p computed from DIN formula

%dsc_list = DSC204_readFiles(['/home/argo/masterarbeit/', ...
%    'DSC204_F1_Phoenix_Messungen/Messungen/Messungen/', ...
%    'ExpDat_16-407-3_mitKorr_*Kmin_H.csv']);

dsc_list = {'ExpDat_16-407-3_mitKorr_20Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_10Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_5Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_2,5Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_1,25Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_0,6Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_0,3Kmin_H.csv'};
            
        
fig1 = figure(1); hold on
ax1 = gca();
fig2 = figure(2); hold on
ax2 = gca();


for i=1:length(dsc_list)

    dsc = DSC204_readFile(dsc_list{i});
    
    q_meas = dsc.data(:,3) ./ dsc.data(:,4) * dsc.mass;
    legend_str = [num2str(dsc.Tinfo.Tstep), ' K/min'];
    plot(ax1, dsc.data(:,1), q_meas, 'DisplayName', legend_str, ...
        'LineWidth', 2.2);
    
    c_p = calc_cp(dsc);
    legend_str = [num2str(dsc.Tinfo.Tstep), ' K/min'];
    plot(ax2, c_p(:,1), c_p(:,2), 'DisplayName', legend_str, ...
        'LineWidth', 1.3);
    
    
    
end

save_root_dir = '/home/argo/masterarbeit/vortrag/images/';

set(fig1, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
set(ax1,'FontSize',25);
xlabel(ax1, 'T_{ref} [°C]');
% ylabel(ax1, '\eta^{\Phi}[mW]');
ylabel(ax1, '\Phi_q^{\eta} [mW]');
set(ax1, 'xlim', [30 160]);
legend(ax1, 'show', 'location', 'northwest');
xlim(ax1, [80, 160]);


print(fig1, [save_root_dir, 'heat_flux_measurement'], '-dpng', '-r200');
close(fig1);


set(fig2, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
set(ax2,'FontSize',25);
xlabel(ax2, 'T [°C]');
ylabel(ax2, 'c_p [mJ/(mg*K)]');
title(ax2, 'DIN 11357 formula')
set(ax2, 'xlim', [30 160]);
legend(ax2, 'show', 'location', 'northwest');
xlim(ax2, [105, 160]);

print(fig2, [save_root_dir, 'c_p_DIN_formula'], '-dpng', '-r200');
close(fig2);


return



%% Plot optimization progress

fit_dir = '2017-12-20_14:25:10_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS_used';

% fit_dir = '2017-12-20_02:01:34_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS';
% fit_dir = '2017-12-17_22:56:22_407_L1=40_L3=0,1_N1=200_N3=50_GN_Gaussians';
% fit_dir = '2017-12-15_00:21:25_407_L1=40_L3=0,1_N1=200_N3=50_GN_FS';
% fit_dir = '2017-12-15_13:42:19_407_L1=40_L3=0,1_N1=200_N3=50_GN_FS';

fit_path = strcat('/home/argo/masterarbeit/fits_data/', fit_dir, '/');

file_list = dir(fit_path);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];


fig = figure(1); clf;
set(fig, 'Units', 'normalized', 'OuterPosition', [0., 0., 0.53, 1.1]); 

ax1 = gca; set(ax1, 'YScale', 'log'); hold on


for i=1:length(nameSubDirs)
        
    clf;
    ax1 = gca; set(ax1, 'YScale', 'log'); hold on

    disp(nameSubDirs{i})
    
    filepath = strcat(fit_path, nameSubDirs{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    num_iterations = length(fit_data.optimization.progress_dx_norm);
    
    plot(ax1, 0:num_iterations, fit_data.optimization.progress_F1_norm, ...
        '-s', 'DisplayName', '||F_1^{(k)}||_2', 'Linewidth', 2.);
    plot(ax1, 0:num_iterations-1, fit_data.optimization.progress_dx_norm, ...
        '-s', 'DisplayName', '||\Deltax^{(k)}||_2', 'Linewidth', 2.);
    plot(ax1, 0:num_iterations-1, fit_data.optimization.progress_t_k, ...
        '-s', 'DisplayName', 't^{(k)}', 'Linewidth', 2.);
    plot(ax1, 0:num_iterations, fit_data.optimization.progress_NOC1, ...
        '-s', 'DisplayName', '||\nablaL^{(k)}||_2', 'Linewidth', 2.);
    
    
    xlim(ax1, [0, num_iterations]);
    xlabel('#Iteration');
    
    heat_rate_str = num2str(fit_data.measurement.dsc_data.Tinfo.Tstep);
    title(ax1, sprintf('Heat rate: %s K/min',heat_rate_str));
    
    legend(ax1, 'show', 'location', 'northeast', 'Orientation', 'horizontal');
    
%     if (i == 4 || i == 6 || i == 7 )
%         legend(ax1, 'show', 'location', 'southwest', 'Orientation', 'vertical');
%     else
%         legend(ax1, 'show', 'location', 'northeast');
%     end
    
    set(ax1,'FontSize',20, 'FontWeight', 'bold');
    grid(ax1, 'on');    
    
    print(fig, [fit_path, nameSubDirs{i}, '/optimization_progress_vortrag'], '-dpng', '-r300');
 
end

