function [enthalpies] = get_fit_results(fit_dir)
%


fit_path = strcat('/home/argo/masterarbeit/fits_data/', fit_dir, '/');

file_list = dir(fit_path);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];


figure();
hold on;
T_domain = 30:0.01:170;
enthalpies = zeros(1,7);
for i=1:length(nameSubDirs)
   
    disp(nameSubDirs{i})
    
    filepath = strcat(fit_path, nameSubDirs{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
    plot(T_domain, c_p, 'DisplayName', strcat(num2str(fit_data.measurement.dsc_data.Tinfo.Tstep), ' K/min'))
    xlabel('T [degC]')
    ylabel('c_p [mJ/(mg K)]')
    
    enthalpies(i) = integral(@(x)c_p_gauss_linear_comb(x, fit_data.optimization.p_optim_end), 100, 160);
    
end
legend('show', 'location', 'northwest');






return



