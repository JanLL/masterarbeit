% comparison of c_p curves for fits with different heat_rate

filename_root = '/home/argo/masterarbeit/fits_data/';
filename_list = ...
    {'2017-08-19_00:28:47_407_20Kmin_lsqnonlin', ...
     '2017-08-17_12:52:24_407_10Kmin_lsqnonlin', ...
     '2017-08-18_17:16:11_407_5Kmin_lsqnonlin', ...
     '2017-08-18_17:38:20_407_2,5Kmin_lsqnonlin', ...
     '2017-08-18_17:52:10_407_1,25Kmin_lsqnonlin', ...
     '2017-08-18_17:57:47_407_0,6Kmin_lsqnonlin', ...
     '2017-08-18_18:06:54_407_0,3Kmin_lsqnonlin'};

num_files = length(filename_list);
 
figure(1);
for n=1:num_files
 
fit_data_path = strjoin(strcat(filename_root, filename_list(n), '/fit_data.mat'));

 
 fit_data = load(fit_data_path);
 p_sim = fit_data.simulation;
 
 T = 30:0.01:160;
 plot(T, p_sim.eval_c_p(T), 'DisplayName', sprintf('heat rate %2.4g K/min', p_sim.heat_rate)); hold on
 
end

legend('show', 'location', 'northwest');
title('Fit results for different heat rates');
xlabel('Temp [degC]');
ylabel('c_p [mJ/(mg K)]')