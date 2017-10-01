%% Position Peak bei variierendem L1

heat_rate = [0.3, 0.6, 1.25, 2.5, 5.0, 10.0];

peak_pos_15 = [129.8, 128.7, 127.2, 125.5, 123.6, 122.0];
peak_pos_12 = [130.0, 129.2, 128.0, 126.5, 124.9, 123.5];
peak_pos_18 = [129.4, 128.3, 126.6, 124.6, 122.5, 120.9];
peak_pos_5  = [131.0, 130.4, 130.0, 129.4, 129.0, 132.2];


figure()
hold on;
plot(heat_rate, peak_pos_15, 'DisplayName', 'L1=15');
plot(heat_rate, peak_pos_12, 'DisplayName', 'L1=12');
plot(heat_rate, peak_pos_18, 'DisplayName', 'L1=18');
plot(heat_rate, peak_pos_5, 'DisplayName', 'L1=5');

legend('show', 'location', 'northoutside')


%% Plot c_p curves for setting L1=15, L3=0.5 for all heat rates
fit_list = {'2017-09-30_15:06:16_407_0.3Kmin_L1=15_L3=0,5', ...
            '2017-09-30_14:55:17_407_0.6Kmin_L1=15_L3=0,5', ...
            '2017-09-30_14:45:27_407_1.25Kmin_L1=15_L3=0,5', ...
            '2017-09-30_14:29:50_407_2.5Kmin_L1=15_L3=0,5', ...
            '2017-09-30_14:19:48_407_5Kmin_L1=15_L3=0,5', ...
            '2017-09-30_14:13:33_407_10Kmin_L1=15_L3=0,5', ...
            '2017-09-30_14:04:55_407_20Kmin_L1=15_L3=0,5'};
        
        
figure();
hold on;
T_domain = 30:0.01:170;
for i=1:length(fit_list)
   
    disp(fit_list{i})
    
    filepath = strcat('/home/argo/masterarbeit/fits_data/', fit_list{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
    plot(T_domain, c_p, 'DisplayName', num2str(fit_data.measurement.dsc_data.Tinfo.Tstep))
    
end
legend('show', 'location', 'northwest');
            



%% Plot c_p curves for setting L1=12, L3=0.5 for all heat rates
fit_list = {'2017-09-29_17:29:00_407_0.3Kmin_L1=12_L3=0,5', ...
            '2017-09-29_17:22:05_407_0.6Kmin_L1=12_L3=0,5', ...
            '2017-09-29_17:14:51_407_1.25Kmin_L1=12_L3=0,5', ...
            '2017-09-29_16:46:36_407_2.5Kmin_L1=12_L3=0,5', ...
            '2017-09-29_16:39:26_407_5Kmin_L1=12_L3=0,5', ...
            '2017-09-29_16:29:22_407_10Kmin_L1=12_L3=0,5'};
        
        
figure();
hold on;
T_domain = 30:0.01:170;
for i=1:length(fit_list)
   
    disp(fit_list{i})
    
    filepath = strcat('/home/argo/masterarbeit/fits_data/', fit_list{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
    plot(T_domain, c_p, 'DisplayName', num2str(fit_data.measurement.dsc_data.Tinfo.Tstep))
    
end
legend('show', 'location', 'northwest');
            



%% Plot c_p curves for setting L1=5, L3=0.5 for all heat rates
fit_list = {'2017-09-30_17:24:48_407_0,3Kmin_L1=5_L3=0,5', ...
            '2017-09-30_17:20:50_407_0,6Kmin_L1=5_L3=0,5', ...
            '2017-09-30_17:13:30_407_1,25Kmin_L1=5_L3=0,5', ...
            '2017-09-30_17:02:45_407_2,5Kmin_L1=5_L3=0,5', ...
            '2017-09-30_16:57:16_407_5Kmin_L1=5_L3=0,5', ...
            '2017-09-30_16:49:07_407_10Kmin_L1=5_L3=0,5', ...
            '2017-09-30_16:38:13_407_20Kmin_L1=5_L3=0,5'};
        
        
figure();
hold on;
T_domain = 30:0.01:170;
for i=1:length(fit_list)
   
    disp(fit_list{i})
    
    filepath = strcat('/home/argo/masterarbeit/fits_data/', fit_list{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
    plot(T_domain, c_p, 'DisplayName', num2str(fit_data.measurement.dsc_data.Tinfo.Tstep))
    
end
legend('show', 'location', 'northwest');




%% Plot c_p curves for setting L1=5, L3=0.3 for all heat rates
fit_list = {'2017-09-30_20:12:36_407_0,3Kmin_L1=5_L3=0,3', ...
            '2017-09-30_20:06:05_407_0,6Kmin_L1=5_L3=0,3', ...
            '2017-09-30_20:00:33_407_1,25Kmin_L1=5_L3=0,3', ...
            '2017-09-30_19:51:57_407_2,5Kmin_L1=5_L3=0,3', ...
            '2017-09-30_19:40:33_407_5Kmin_L1=5_L3=0,3', ...
            '2017-09-30_19:28:48_407_10Kmin_L1=5_L3=0,3', ...
            '2017-09-30_19:14:21_407_20Kmin_L1=5_L3=0,3'};
        
        
figure();
hold on;
T_domain = 30:0.01:170;
for i=1:length(fit_list)
   
    disp(fit_list{i})
    
    filepath = strcat('/home/argo/masterarbeit/fits_data/', fit_list{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
    plot(T_domain, c_p, 'DisplayName', num2str(fit_data.measurement.dsc_data.Tinfo.Tstep))
    
end
legend('show', 'location', 'northwest');
            
            
            
%% Plot c_p curves for setting L1=5, L3=0.7 for all heat rates
fit_list = {'2017-09-30_21:19:55_407_0,3Kmin_L1=5_L3=0,7', ...
            '2017-09-30_21:10:03_407_0,6Kmin_L1=5_L3=0,7', ...
            '2017-09-30_21:04:16_407_1,25Kmin_L1=5_L3=0,7', ...
            '2017-09-30_20:56:24_407_2,5Kmin_L1=5_L3=0,7', ...
            '2017-09-30_20:51:57_407_5Kmin_L1=5_L3=0,7', ...
            '2017-09-30_20:43:43_407_10Kmin_L1=5_L3=0,7', ...
            '2017-09-30_20:38:09_407_20Kmin_L1=5_L3=0,7'};
        
figure();
hold on;
T_domain = 30:0.01:170;
for i=1:length(fit_list)
   
    disp(fit_list{i})
    
    filepath = strcat('/home/argo/masterarbeit/fits_data/', fit_list{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
    plot(T_domain, c_p, 'DisplayName', num2str(fit_data.measurement.dsc_data.Tinfo.Tstep))
    
end
legend('show', 'location', 'northwest');
            
            
            


%% Plot c_p curves for setting L1=5, L3=0.1 for all heat rates
fit_list = {'2017-10-01_11:03:29_407_0,3Kmin_L1=5_L3=0,1', ...
            '2017-10-01_11:00:22_407_0,6Kmin_L1=5_L3=0,1', ...
            '2017-10-01_10:55:12_407_1,25Kmin_L1=5_L3=0,1', ...
            '2017-10-01_10:49:35_407_2,5Kmin_L1=5_L3=0,1', ...
            '2017-10-01_10:42:55_407_5Kmin_L1=5_L3=0,1', ...
            '2017-10-01_10:31:15_407_10Kmin_L1=5_L3=0,1', ...
            '2017-10-01_10:17:21_407_20Kmin_L1=5_L3=0,1'};
        
figure();
hold on;
T_domain = 30:0.01:170;
for i=1:length(fit_list)
   
    disp(fit_list{i})
    
    filepath = strcat('/home/argo/masterarbeit/fits_data/', fit_list{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
    plot(T_domain, c_p, 'DisplayName', num2str(fit_data.measurement.dsc_data.Tinfo.Tstep))
    
end
legend('show', 'location', 'northwest');
            


%% Plot c_p curves for setting L1=15, L3=0.5 and SILVER for all heat rates
fit_list = {'2017-10-01_13:09:55_407_0,3Kmin_L1=15_L3=0,5', ...
            '2017-10-01_13:05:45_407_0,6Kmin_L1=15_L3=0,5', ...
            '2017-10-01_12:50:60_407_2,5Kmin_L1=15_L3=0,5', ...
            '2017-10-01_12:44:35_407_5Kmin_L1=15_L3=0,5', ...
            '2017-10-01_12:30:28_407_10Kmin_L1=15_L3=0,5', ...
            '2017-10-01_12:20:23_407_20Kmin_L1=15_L3=0,5'};
        
figure();
hold on;
T_domain = 30:0.01:170;
for i=1:length(fit_list)
   
    disp(fit_list{i})
    
    filepath = strcat('/home/argo/masterarbeit/fits_data/', fit_list{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
    plot(T_domain, c_p, 'DisplayName', num2str(fit_data.measurement.dsc_data.Tinfo.Tstep))
    
end
legend('show', 'location', 'northwest');
            
            
