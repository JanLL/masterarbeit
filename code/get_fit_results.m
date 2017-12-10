function [fig1, enthalpies] = get_fit_results(fit_dir)
%


fit_path = strcat('/home/argo/masterarbeit/fits_data/', fit_dir, '/');

file_list = dir(fit_path);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];

% Figure for all c_p plots
fig1 = figure(66); clf; ax1 = gca;
hold on;

% Figure for integration process plots
fig2 = figure(67); clf;
ax2 = {subplot(4,2,1), subplot(4,2,2), subplot(4,2,3), subplot(4,2,4), ...
       subplot(4,2,5), subplot(4,2,6), subplot(4,2,7)};
for i=1:7
    set(ax2{i}, 'YScale', 'log');
    hold(ax2{i});
end
   
T_domain = 30:0.01:170;
enthalpies = zeros(1,7);

for i=1:length(nameSubDirs)
    
    disp(nameSubDirs{i})
    
    filepath = strcat(fit_path, nameSubDirs{i}, '/fit_data.mat');
    fit_data = load(filepath);
            
    % c_p plot
    if (strcmp(fit_data.optimization.c_p_param_type, 'gauss_linear_comb'))
        c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
        enthalpies(i) = integral(@(x)c_p_gauss_linear_comb(x, fit_data.optimization.p_optim_end), 100, 160);
    elseif (strcmp(fit_data.optimization.c_p_param_type, 'fraser_suzuki'))
        c_p = c_p_fs(T_domain, fit_data.optimization.p_optim_end);
        enthalpies(i) = integral(@(x)c_p_fs(x, fit_data.optimization.p_optim_end), 110, 160);
    end
    
    plot(ax1, T_domain, c_p, 'DisplayName', strcat(num2str(fit_data.measurement.dsc_data.Tinfo.Tstep), ' K/min'))
    xlabel('T [degC]')
    ylabel('c_p [mJ/(mg K)]')
    
    % optimization progress plot
    %plot(ax2{i}, fit_data.optimization.progress_F1_norm, 'DisplayName', '||F1||_2');
%     plot(ax2{i}, fit_data.optimization.progress_dx_norm, 'DisplayName', '||dx||_2');
%     plot(ax2{i}, fit_data.optimization.progress_t_k, 'DisplayName', 't_k');
    
end

legend(ax1, 'show', 'location', 'northwest');





return



