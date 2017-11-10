function [fig, enthalpies] = get_fit_results(fit_dir)
%


fit_path = strcat('/home/argo/masterarbeit/fits_data/', fit_dir, '/');

file_list = dir(fit_path);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];


fig = figure();
hold on;
T_domain = 30:0.01:170;
enthalpies = zeros(1,7);

for i=1:length(nameSubDirs)
   
    disp(nameSubDirs{i})
    
    filepath = strcat(fit_path, nameSubDirs{i}, '/fit_data.mat');
    fit_data = load(filepath);
        
    if (strcmp(fit_data.optimization.c_p_param_type, 'gauss_linear_comb'))
        c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
        enthalpies(i) = integral(@(x)c_p_gauss_linear_comb(x, fit_data.optimization.p_optim_end), 100, 160);
    elseif (strcmp(fit_data.optimization.c_p_param_type, 'fraser_suzuki'))
        c_p = c_p_fs(T_domain, fit_data.optimization.p_optim_end);
        enthalpies(i) = integral(@(x)c_p_fs(x, fit_data.optimization.p_optim_end), 110, 160);
    end
    
    plot(T_domain, c_p, 'DisplayName', strcat(num2str(fit_data.measurement.dsc_data.Tinfo.Tstep), ' K/min'))
    xlabel('T [degC]')
    ylabel('c_p [mJ/(mg K)]')
    
    %
    
end
legend('show', 'location', 'northwest');






return



