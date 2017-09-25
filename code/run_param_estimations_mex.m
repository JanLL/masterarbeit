
% Simulation
simulation = struct();
simulation.L1 = 15.;
simulation.L3 = 0.5;
simulation.N1 = 200;
simulation.N3 = 50;
 
simulation.lambda_Const = 23.; 
simulation.rho_Const = 8.9;    
simulation.c_p_Const = 0.41;   

simulation.lambda_pcm = 0.96;  
simulation.rho_pcm = 0.85;    

simulation.T_0 = 10.;
simulation.T_end = 200.;

% Measurements
dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

% Optimization
optimization = struct();
optimization.c_p_param_type = 'gauss_linear_comb';
optimization.start_values = [37.2793,    1.2498,  123.0044, ...
                             1.1022,   13.2169,  104.8824, ...
                             -0.3,    3.7355,  124.3034, ...
                             2.8868, 7.1009,  118.2194, ...
                             0.5846,   30.5770,   76.8653, ...
                             0.1, 1., 125., ...
                             0.1, 1., 127., ...
                             0.1, 1., 129., ...
                             0.1, 1., 70., ...
                             0.1, 1., 80., ...
                             1., 1.6483];
num_opt_params = length(optimization.start_values);

optimization.p_optim_estimable = true(length(optimization.start_values), 1);

optimization.lb = ones(num_opt_params, 1) * -1.;
optimization.ub = ones(num_opt_params, 1) * 500.;


% General options
options.init_value_test = false;

% Run optimization #1
fit_data = param_estimation_mex2(simulation, dsc_measurement, optimization, options);







