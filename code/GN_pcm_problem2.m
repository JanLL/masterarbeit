function [fit_data] = GN_pcm_problem2(simulation, dsc_measurement, optimization, options)
% Description!!!

% Simulation parameters
L1 = simulation.L1;
L3 = simulation.L3;
N1 = simulation.N1;
N3 = simulation.N3; 

 
lambda_Const = simulation.lambda_Const;
rho_Const    = simulation.rho_Const;
c_p_Const    = simulation.c_p_Const;

lambda_pcm = simulation.lambda_pcm;

T_0       = simulation.T_0;
T_end     = simulation.T_end;


% Grid generation
n_tr = simulation.grid_n_tr;
n_m = simulation.grid_n_m;
t = simulation.grid_t;

spatial_gridsize = build_grid(N1, N3, L1, L3, n_tr, n_m, t);


% Measurement
m_pcm = dsc_measurement.mass;
heat_rate = dsc_measurement.Tinfo.Tstep;

switch heat_rate
    case 20
        heat_rate = 19.97;
    case 10
        heat_rate = 9.98;
    case 5
        heat_rate = 4.99;
    case 2.5
        heat_rate = 2.496;
    case 1.25
        heat_rate = 1.2476;
    case 0.6
        heat_rate = 0.5989;
    case 0.3
        heat_rate = 0.2994;
    otherwise
        error('Heat rate modification failed!')
end

simulation.heat_rate = heat_rate;


% TODO: sinnvolles Intervall automatisch waehlen ... wobei das hier fuer
% alle Messungen bisher ganz gut war
index_T_dsc = [find(dsc_measurement.data(:,1) > 29, 1, 'first'), ...
                           find(dsc_measurement.data(:,1) < 157.9, 1, 'last')];

T_ref_dsc = dsc_measurement.data(index_T_dsc(1):index_T_dsc(2),1);
q_dsc = (dsc_measurement.data(index_T_dsc(1):index_T_dsc(2),3) ...
      ./ dsc_measurement.data(index_T_dsc(1):index_T_dsc(2),4)) * m_pcm;

num_meas = length(q_dsc);


% Some pre calculations
N = N1+N3;
a_Const = lambda_Const/(rho_Const*c_p_Const);
heat_rate_s = heat_rate / 60; % [K/min] -> [K/s]

% Analytical solution for T_ref: Solve non-linear system to get time t 
% where T_ref (temp. at crucible) reaches T_ref_meas (time where heat flux
% measurement was done)
n = 100;
meas_data = zeros(num_meas, 2);
for i=1:length(T_ref_dsc(:,1))

    F = @(t) analytical_sol(L1,t,n,T_0, heat_rate_s, a_Const) - T_ref_dsc(i,1);

    t_guess = (T_ref_dsc(i,1) - T_0)/heat_rate_s;
    
    fsolve_options = optimoptions('fsolve','Display','none');
    meas_data(i,1) = fsolve(F, t_guess, fsolve_options);
end
meas_data(:,2) = q_dsc;

% Optimization
p_optim_fixed = optimization.start_values(~optimization.p_optim_estimable);

% Initialize Mex File
sim_params_vec = [L1, L3, N1, N3, lambda_Const, rho_Const, c_p_Const, ...
                  lambda_pcm, m_pcm, ...
                  heat_rate, T_0, T_end];

heat1D_pcm('reset');
heat1D_pcm('init', sim_params_vec, spatial_gridsize, meas_data, optimization.c_p_param_type);


% Initialize figures for plots during optimization
figure(1); % q_pcm_in plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();
figure(3); % dqdp plot
ax3 = gca();

F1_func = @(p_optim) compute_q_dqdp_mex(...
    p_optim, optimization.p_optim_estimable, optimization.p_optim_fixed, ...
    optimization.c_p_param_type, T_ref_dsc, q_dsc, ax1, ax2, ax3);
F2_func = @(p) GN_test_fct_F2(p);


%%%%%%%%%%%%%% INITIAL VALUE TEST %%%%%%%%%%%%%%%%%%%%
if options.init_value_test
    [res, Jac] = F1_func(optimization.start_values(optimization.p_optim_estimable));
    fit_data = false;
  
    return;
end


% Start optimizing
[p_optim_end, optim_output] = ...
    GN_ass(F1_func, ...
           F2_func, ...
           optimization.start_values(optimization.p_optim_estimable), ...
           optimization.lb(optimization.p_optim_estimable), ...
           optimization.ub(optimization.p_optim_estimable), ...
           options.GN_options);
       
p_optim_all = zeros(1,length(optimization.p_optim_estimable));
p_optim_all(optimization.p_optim_estimable) = p_optim_end;
p_optim_all(~optimization.p_optim_estimable) = p_optim_fixed;

%%%%%%%%%%%%%%%%% SAVE FIT RESULTS %%%%%%%%%%%%%%%%%%
save_path_root = options.save_path_root;

% optim_progress = [optim_output.progress_F1_norm, ...
%                   optim_output.progress_dx_norm, ...
%                   optim_output.progress_t_k];

fit_data = GN_save_fit(...
    save_path_root, simulation, dsc_measurement, index_T_dsc, ...
    optimization, p_optim_all, optim_output);


end

