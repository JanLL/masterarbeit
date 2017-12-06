% Open error log file
errLogFileID = fopen('/home/argo/masterarbeit/ErrLog.txt', 'a');
fprintf(errLogFileID, '\n\nStart new Set of parameter estimations at %s...\n', datetime('now'));


% Simulation
simulation = struct();
simulation.L1 = 40.;
simulation.L3 = 0.1;
simulation.N1 = 300;
simulation.N3 = 50; 

% Silver
simulation.lambda_Const = 430.; 
simulation.rho_Const = 10.49;    
simulation.c_p_Const = 0.235; 

simulation.lambda_pcm = 0.96;  

simulation.T_0 = 10.;
simulation.T_end = 200.;

simulation.grid_n_tr = 0.1;
simulation.grid_n_m = 0.01;
simulation.grid_t = 0.999;


% Measurements
dsc_filename = 'ExpDat_16-407-3_mitKorr_20Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Optimization
optimization = struct();

% optimization.c_p_param_type = 'fraser_suzuki';
optimization.c_p_param_type = 'gauss_linear_comb';

%%%%%%%%%% Set optimization variables Fraser Suzuki %%%%%%%%%%%%%%%%%%%%%%
if (strcmp(optimization.c_p_param_type, 'fraser_suzuki'))
    
    c_p_type = 'FS';
    
    h  =  10.0;
    r  =  2.0;
    wr =  5.0;
    sr =   0.3;
    z  = 130.0;
    m  = 0.1;
    b  =   3.0;
    p_fraser_suzuki = [h, r, wr, sr, z, m, b].';

    optimization.start_values = p_fraser_suzuki;

    % choose free(true)/fixed(false) parameters to optimize
    optimization.p_optim_estimable = true(length(optimization.start_values), 1);
    optimization.p_optim_estimable(2) = false;  % fix "r"
    optimization.p_optim_fixed = optimization.p_optim_start(~optimization.p_optim_estimable);

    optimization.lb = -inf*ones(1,length(optimization.start_values));
    optimization.ub = +inf*ones(1,length(optimization.start_values));

elseif (strcmp(optimization.c_p_param_type, 'gauss_linear_comb'))
    
    c_p_type = 'Gaussians';
    
%     fit_data = load(['/home/argo/masterarbeit/fits_data/', ...
%                      '2017-12-03_19:49:07_407_L1=40_L3=0.1_N1=500_N3=50/', ...
%                      '2017-12-03_20:09:10_407_20Kmin_L1=40_L3=0,1/fit_data.mat']);
%     optimization.p_optim_start = fit_data.optimization.p_optim_end(1:32).';
    
    % 5 Gausse aktiviert
    optimization.start_values = ones(32,1);
    optimization.start_values(1:3:15) = [0.5, 0.5, 0.5, 20., -2.];
    optimization.start_values(2:3:15) = [2., 2., 2., 5., 4.];
    optimization.start_values(3:3:15) = [110., 115., 120., 128., 145.];
    
    optimization.start_values(31) = 1.07;  % linear
    optimization.start_values(32) = 1.78;  % const     

    % choose free(true)/fixed(false) parameters to optimize
    optimization.p_optim_estimable = true(length(optimization.start_values), 1);
    optimization.p_optim_estimable(16:30) = false;  % deactivate last 5 Gaussians
    
    
    optimization.p_optim_fixed = optimization.start_values(~optimization.p_optim_estimable);
    
    optimization.lb = -inf*ones(1,length(optimization.start_values));
    optimization.lb(2:3:30) = 0.5;
    
    optimization.ub = +inf*ones(1,length(optimization.start_values));
    optimization.ub(2:3:30) = 100.;
    
end


% General options
options = struct();
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
                            'N3=', N3_str, '_', ...
                            'GN_', c_p_type, ...
                            '/');

if exist(save_path_root_str, 'dir') == 0
    disp('Creating root directory');
    mkdir(save_path_root_str);
end

options.save_path_root = save_path_root_str;

GN_options = struct;
GN_options.decomposition = 'SVD';
GN_options.TOL_ineq = 1e-8;  % constraint active when -TOL < F_3i < TOL

% Termination criteria
GN_options.TOL_dx_norm = 1e-8;
GN_options.TOL_t_k = 1e-7;
GN_options.max_iterations = 5;
GN_options.t_k_start = 0.3;

options.GN_options = GN_options;

% Run optimization #1
try
    fit_data = GN_pcm_problem2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
    fprintf('%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end


% Optimization #2 settings
optimization.start_values = fit_data.optimization.p_optim_end;
dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Run optimization #2
try
    fit_data = GN_pcm_problem2(simulation, dsc_measurement, optimization, options);
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
    fit_data = GN_pcm_problem2(simulation, dsc_measurement, optimization, options);
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
    fit_data = GN_pcm_problem2(simulation, dsc_measurement, optimization, options);
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
    fit_data = GN_pcm_problem2(simulation, dsc_measurement, optimization, options);
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
    fit_data = GN_pcm_problem2(simulation, dsc_measurement, optimization, options);
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
    fit_data = GN_pcm_problem2(simulation, dsc_measurement, optimization, options);
catch Err
    fprintf(errLogFileID, '%s\tError %s occured with heat_rate=%2.4g.\n', ...
        datetime('now'), Err.identifier, dsc_measurement.Tinfo.Tstep);
end



fprintf(errLogFileID, 'Finished Set of parameter estimations at %s...\n', datetime('now'));
fclose(errLogFileID);
return;



