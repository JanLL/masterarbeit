% Open error log file
errLogFileID = fopen('/home/argo/masterarbeit/ErrLog.txt', 'a');
fprintf(errLogFileID, '\n\nStart new Set of parameter estimations at %s...\n', datetime('now'));


% Simulation
simulation = struct();
simulation.L1 = 15.;
simulation.L3 = 0.5;
simulation.N1 = 200;
simulation.N3 = 50;

% Constantan
simulation.lambda_Const = 23.; 
simulation.rho_Const = 8.9;    
simulation.c_p_Const = 0.41;  

% Silver
% simulation.lambda_Const = 430.; 
% simulation.rho_Const = 10.49;    
% simulation.c_p_Const = 0.235; 

simulation.lambda_pcm = 0.96;  

simulation.T_0 = 10.;
simulation.T_end = 200.;

% Measurements
dsc_filename = 'ExpDat_16-407-3_mitKorr_20Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Optimization
optimization = struct();
optimization.c_p_param_type = 'gauss_linear_comb';

fit_data = load('/home/argo/masterarbeit/fits_data/2017-09-30_16:49:07_407_10Kmin_L1=5_L3=0,5/fit_data.mat');
optimization.start_values = fit_data.optimization.p_optim_end;

% optimization.start_values = [10.,    1.,  123.0044, ...
%                              0.,    1.,  104.8824, ...
%                              0.,    1.,  124.3034, ...
%                              0.,    1.,  118.2194, ...
%                              0.,    1.,   76.8653, ...
%                              0.1, 1., 125., ...
%                              0.1, 1., 127., ...
%                              0.1, 1., 129., ...
%                              -1., 1., 133., ...
%                              0.1, 1., 80., ...
%                              0.5, 3.];
num_opt_params = length(optimization.start_values);

optimization.p_optim_estimable = true(length(optimization.start_values), 1);

optimization.lb = zeros(num_opt_params,1);
optimization.lb(1:3:30) = -2.;
optimization.lb(2:3:30) = 0.15;

optimization.ub = ones(num_opt_params,1) * inf;


% General options
options.init_value_test = false;


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



