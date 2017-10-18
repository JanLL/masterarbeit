% Open error log file
errLogFileID = fopen('/home/argo/masterarbeit/ErrLog.txt', 'a');
fprintf(errLogFileID, '\n\nStart new Set of parameter estimations at %s...\n', datetime('now'));


% Simulation
simulation = struct();
simulation.L1 = 40.;
simulation.L3 = 0.1;
simulation.N1 = 500;
simulation.N3 = 50;

% Constantan
% simulation.lambda_Const = 23.; 
% simulation.rho_Const = 8.9;    
% simulation.c_p_Const = 0.41;  

% Silver
simulation.lambda_Const = 430.; 
simulation.rho_Const = 10.49;    
simulation.c_p_Const = 0.235; 

simulation.lambda_pcm = 0.96;  

simulation.T_0 = 10.;
simulation.T_end = 200.;

% Measurements
dsc_filename = 'ExpDat_16-407-3_mitKorr_20Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Optimization
optimization = struct();
optimization.c_p_param_type = 'gauss_linear_comb';

fit_data = load(['/home/argo/masterarbeit/fits_data/', ...
                 '2017-10-16_04:45:00_407_L1=40_L3=0.1_N1=200_N3=50/', ...
                 '2017-10-16_04:45:34_407_20Kmin_L1=40_L3=0,1/fit_data.mat']);

optimization.start_values = fit_data.optimization.p_optim_end;
% optimization.start_values = [10.,   1.,  123.0, ...
%                              1.,    1.,  115.0, ...
%                              1.,    1.,  125.0, ...
%                              1.,    1.,  127.0, ...
%                              -0.1,    1.,  129.0, ...
%                              0., 1., 0., ...
%                              0., 1., 0., ...
%                              0., 1., 0., ...
%                              0., 1., 0., ...
%                              0., 1., 0., ...
%                              0.1, 2.];
num_opt_params = length(optimization.start_values);

optimization.p_optim_estimable = true(length(optimization.start_values), 1);
%optimization.p_optim_estimable(16:30) = false;  %  fix 5 Gaussians

optimization.lb = zeros(num_opt_params,1);
optimization.lb(1:3:30) = -2.;
optimization.lb(2:3:30) = 0.15;

optimization.ub = ones(num_opt_params,1) * inf;


% General options
options.init_value_test = false;

datetime_cell = num2cell(int32(clock));
datetime_str = sprintf('%04i-%02i-%02i_%02i:%02i:%02i', datetime_cell{:});

mass_code_str = dsc_measurement.fileSpec(11:13);
L1_str = num2str(simulation.L1);
L3_str = num2str(simulation.L3);
N1_str = num2str(simulation.N1);
N3_str = num2str(simulation.N3);

save_path_root_str = strcat('/home/argo/masterarbeit/fits_data/', ...
                            datetime_str, '_', ...
                            mass_code_str, '_', ...
                            'L1=', L1_str, '_', ...
                            'L3=', L3_str, '_', ...
                            'N1=', N1_str, '_', ...
                            'N3=', N3_str, ...
                            '/');

if exist(save_path_root_str, 'dir') == 0
    disp('Creating root directory');
    mkdir(save_path_root_str);
end

options.save_path_root = save_path_root_str;


% Run optimization #1
try
    fit_data = param_estimation_mex2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end




% Optimization #2 settings
optimization.start_values = fit_data.optimization.p_optim_end;
dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Run optimization #2
try
    fit_data = param_estimation_mex2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end




% Optimization #3 settings
optimization.start_values = fit_data.optimization.p_optim_end;
dsc_filename = 'ExpDat_16-407-3_mitKorr_5Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Run optimization #3
try
    fit_data = param_estimation_mex2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end


% Optimization #4 settings
optimization.start_values = fit_data.optimization.p_optim_end;
dsc_filename = 'ExpDat_16-407-3_mitKorr_2,5Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Run optimization #4
try
    fit_data = param_estimation_mex2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end


% Optimization #5 settings
optimization.start_values = fit_data.optimization.p_optim_end;
dsc_filename = 'ExpDat_16-407-3_mitKorr_1,25Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Run optimization #5
try
    fit_data = param_estimation_mex2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end


% Optimization #6 settings
optimization.start_values = fit_data.optimization.p_optim_end;
dsc_filename = 'ExpDat_16-407-3_mitKorr_0,6Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Run optimization #6
try
    fit_data = param_estimation_mex2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end


% Optimization #7 settings
optimization.start_values = fit_data.optimization.p_optim_end;
dsc_filename = 'ExpDat_16-407-3_mitKorr_0,3Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Run optimization #7
try
    fit_data = param_estimation_mex2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end



fprintf(errLogFileID, 'Finished Set of parameter estimations at %s...\n', datetime('now'));
fclose(errLogFileID);
return;



