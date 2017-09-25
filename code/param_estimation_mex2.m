function [fit_data] = param_estimation_mex2(simulation, dsc_measurement, optimization, options)
% 



% Simulation parameters
L1 = simulation.L1;
L3 = simulation.L3;
N1 = simulation.N1;
N3 = simulation.N3; 
 
lambda_Const = simulation.lambda_Const; 
rho_Const    = simulation.rho_Const;    
c_p_Const    = simulation.c_p_Const;   

lambda_pcm = simulation.lambda_pcm;  
rho_pcm    = simulation.rho_pcm;    

T_0       = simulation.T_0;
T_end     = simulation.T_end;
              
% Measurement
m_pcm = dsc_measurement.mass;
heat_rate = dsc_measurement.Tinfo.Tstep;

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
                  lambda_pcm, rho_pcm, m_pcm, ...
                  heat_rate, T_0, T_end];

heat1D_pcm('reset');
heat1D_pcm('init', sim_params_vec, meas_data, optimization.c_p_param_type);


% Initialize figures for plots during optimization
figure(1); % q_pcm_in plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();
figure(3); % dqdp plot
ax3 = gca();

compute_q_dqdp_mex_expl = @(p_optim) compute_q_dqdp_mex(...
    p_optim, optimization.p_optim_estimable, p_optim_fixed, ...
    optimization.c_p_param_type, T_ref_dsc, q_dsc, ax1, ax2, ax3);


%%%%%%%%%%%%%% INITIAL VALUE TEST %%%%%%%%%%%%%%%%%%%%
if options.init_value_test
    [res, Jac] = compute_q_dqdp_mex_expl(optimization.start_values);
    fit_data = false;
    return;
end


%%%%%%%%%%%%%% SOLVE OPTIMIZATION PROBLEM %%%%%%%%%%%%
global p_optim_process;
opt_constraints = {optimization.lb, optimization.ub};

opt_options = optimoptions('lsqnonlin', ...
                           'Display', 'iter-detailed', ...
                           'OutputFcn', @disp_aux, ...
                           'SpecifyObjectiveGradient', true, ...
                           'MaxIter', 5);


tic;
[p_optim,~,~,~,optim_output,~,jac_output] = lsqnonlin(...
    compute_q_dqdp_mex_expl, ...
    optimization.start_values(optimization.p_optim_estimable), ...
    opt_constraints{:}, ...
    opt_options);
optim_duration = toc;
fprintf('Optimization took %2.2f seconds.\n', optim_duration);

%%%%%%%%%%%%%%%%% SAVE FIT RESULTS %%%%%%%%%%%%%%%%%%
save_path_root = '/home/argo/masterarbeit/fits_data/';
[residuum_end, ~] = heat1D_pcm('optimization', p_optim);

fit_data = save_fit_mex(...
    save_path_root, simulation, dsc_measurement, index_T_dsc, ...
    optimization, p_optim, optim_output, optim_duration, ...
    residuum_end, jac_output, p_optim_process);


end

