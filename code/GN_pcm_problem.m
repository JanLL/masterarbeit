% Standalone test for c_p optimization in PCM problem to test Gauss-Newton-
% Active Set Strategy.


%%%%%%%%%% Get measurement data %%%%%%%%%%%%%%%%%%%
dsc_filename = 'ExpDat_16-407-3_mitKorr_20Kmin_H.csv';
dsc = DSC204_readFile(dsc_filename);

m_pcm = dsc.mass;

% TODO: sinnvolles Intervall automatisch waehlen ... wobei das hier fuer
% alle Messungen bisher ganz gut war
index_T_dsc = [find(dsc.data(:,1) > 29, 1, 'first'), ...
               find(dsc.data(:,1) < 157.9, 1, 'last')];

T_ref_dsc = dsc.data(index_T_dsc(1):index_T_dsc(2),1);
q_dsc = (dsc.data(index_T_dsc(1):index_T_dsc(2),3) ...
      ./ dsc.data(index_T_dsc(1):index_T_dsc(2),4)) * m_pcm;

num_meas = length(q_dsc);


%%%%%%%%%% Set Simulation parameters %%%%%%%%%%%%%%%%%
L1 = 40;  % [mm]
L3 = 0.1;  % [mm]

% n_pcm = 0.2;
% N = 350;
% N3 = N*n_pcm;
% N1 = N - N3;

N3 = 50;  % error if N3=0
N1 = 300;
N = N1 + N3;

% Constantan
% lambda_Const = 23.;  % [mW/(mm*K)]
% rho_Const = 8.9;     % [mg/mm^3]
% c_p_Const = 0.41;    % [mJ/(mg*K)]

% Silver
lambda_Const = 430.; 
rho_Const = 10.49;    
c_p_Const = 0.235; 

lambda_pcm = 0.96;   % [mW/(mm*K)]
rho_pcm = 0.85;      % [mg/mm^3]

heat_rate = dsc.Tinfo.Tstep;     % [K/min]
T_0 = 10.;           % Start temperature oven [degC]
T_end = 200.;        % End temperature oven [degC]

sim_params_vec = [L1, L3, N1, N3, lambda_Const, rho_Const, c_p_Const, ...
                  lambda_pcm, m_pcm, ...
                  heat_rate, T_0, T_end];

heat_rate_s = heat_rate / 60; % [K/min] -> [K/s]


% Grid generation NEW
n_pcm = N3 / (N1 + N3);
n_tr = 0.1;
n_m = 0.01;
t = 0.999;

N_pcm = N * n_pcm;
N_tr = N * n_tr;
N_m = N * n_m;

b = N-1 - N_pcm - N_m - N_tr/2;
gamma = 2/N_tr * log(t/(1-t));

W11 = sum(1./(exp(gamma*((0:N1-2) - b)) + 1));
W12 = N1 - 1 - sum(1./(exp(gamma*((0:N1-2) - b)) + 1));
W21 = sum(1./(exp(gamma*((0:N1+N3-2) - b)) + 1));
W22 = N1 + N3 - 1 - sum(1./(exp(gamma*((0:N1+N3-2) - b)) + 1));

W = [W11, W12;
     W21, W22];

 dx = W\[L1;L1+L3];
 dx_Ag = dx(1);
 dx_pcm = dx(2);
 
dchi = @(x_tilde) (dx_Ag - dx_pcm) ./ (exp(gamma*(x_tilde-b)) + 1) + dx_pcm;

spatial_gridsize = dchi(0:N-2)';

% sum(spatial_gridsize(1:N1-1))
% sum(spatial_gridsize(1:N1+N3-1))
% figure(10)
% plot(spatial_gridsize, 'x')
% return




%%%%%%%%%% some pre-calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_Const = lambda_Const/(rho_Const*c_p_Const);

% Analytical solution for T_ref: Solve non-linear system to get time t 
% where T_ref (temp. at crucible) reaches T_ref_meas (time where heat flux
% measurement was done)
n = 100;
a = lambda_Const / (c_p_Const * rho_Const);

meas_data = zeros(num_meas, 2);
for i=1:length(T_ref_dsc(:,1))

    F = @(t) analytical_sol(L1,t,n,T_0, heat_rate_s, a) - T_ref_dsc(i,1);

    t_guess = (T_ref_dsc(i,1) - T_0)/heat_rate_s;
    
    fsolve_options = optimoptions('fsolve','Display','none');
    meas_data(i,1) = fsolve(F, t_guess, fsolve_options);
end

meas_data(:,2) = q_dsc;

% optimization.c_p_param_type = 'fraser_suzuki';
optimization.c_p_param_type = 'gauss_linear_comb';

%%%%%%%%%% Set optimization variables Fraser Suzuki %%%%%%%%%%%%%%%%%%%%%%
if (strcmp(optimization.c_p_param_type, 'fraser_suzuki'))
    
    h  =  10.0;
    r  =  2.0;
    wr =  5.0;
    sr =   0.3;
    z  = 130.0;
    m  = 0.1;
    b  =   3.0;
    p_fraser_suzuki = [h, r, wr, sr, z, m, b].';

    optimization.p_optim_start = p_fraser_suzuki;

    % choose free(true)/fixed(false) parameters to optimize
    optimization.p_optim_estimable = true(length(optimization.p_optim_start), 1);
    optimization.p_optim_estimable(2) = false;  % fix "r"
    optimization.p_optim_fixed = optimization.p_optim_start(~optimization.p_optim_estimable);

    optimization.lb = -inf*ones(1,length(optimization.p_optim_start));
    optimization.ub = +inf*ones(1,length(optimization.p_optim_start));

elseif (strcmp(optimization.c_p_param_type, 'gauss_linear_comb'))
    
%     fit_data = load(['/home/argo/masterarbeit/fits_data/', ...
%                      '2017-12-03_19:49:07_407_L1=40_L3=0.1_N1=500_N3=50/', ...
%                      '2017-12-03_20:09:10_407_20Kmin_L1=40_L3=0,1/fit_data.mat']);
%     optimization.p_optim_start = fit_data.optimization.p_optim_end(1:32).';
    
    % 4 Gausse aktiviert
    optimization.p_optim_start = ones(32,1);
    optimization.p_optim_start(1:3:15) = [0.5, 0.5, 0.5, 20., -2.];
    optimization.p_optim_start(2:3:15) = [2., 2., 2., 5., 4.];
    optimization.p_optim_start(3:3:15) = [110., 115., 120., 128., 145.];
    
    optimization.p_optim_start(31) = 1.07;  % linear
    optimization.p_optim_start(32) = 1.78;  % const  

%     % 10 Gausses
%     optimization.p_optim_start = ones(32,1);
%     optimization.p_optim_start(1:3:30) = ...
%         [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 10., -1., -1., -1.];
%     optimization.p_optim_start(2:3:30) = ...
%         [2., 2., 2., 2., 2., 2., 5., 2., 2., 2.];
%     optimization.p_optim_start(3:3:30) = ...
%         [115, 118, 120, 123, 126, 128, 130, 135, 138, 142];
%      
%     optimization.p_optim_start(31) = 1.07;  % linear
%     optimization.p_optim_start(32) = 1.78;  % const    

    % choose free(true)/fixed(false) parameters to optimize
    optimization.p_optim_estimable = true(length(optimization.p_optim_start), 1);
    optimization.p_optim_estimable(16:30) = false;  % deactivate last 5 Gaussians
    
    
    optimization.p_optim_fixed = optimization.p_optim_start(~optimization.p_optim_estimable);
    
    optimization.lb = -inf*ones(1,length(optimization.p_optim_start));
    optimization.lb(2:3:30) = 0.5;
    
    optimization.ub = +inf*ones(1,length(optimization.p_optim_start));
    optimization.ub(2:3:30) = 100.;
    
end


heat1D_pcm('reset');
heat1D_pcm('init', sim_params_vec, spatial_gridsize, meas_data, optimization.c_p_param_type);


figure(1); % q_pcm_in plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();
figure(3); % dqdp plot
ax3 = gca();

F1_func = @(p_optim) compute_q_dqdp_mex(...
    p_optim, optimization.p_optim_estimable, optimization.p_optim_fixed, optimization.c_p_param_type, T_ref_dsc, q_dsc, ax1, ax2, ax3);
F2_func = @(p) GN_test_fct_F2(p);

%%%%%%%%%%%%%% INITIAL VALUE TEST %%%%%%%%%%%%%%%%%%%%
% [F1, J1] = F1_func(optimization.p_optim_start(optimization.p_optim_estimable));
% return


GN_options = struct;
GN_options.decomposition = 'SVD';
GN_options.TOL_ineq = 1e-8;  % constraint active when -TOL < F_3i < TOL

% Termination criteria
GN_options.TOL_dx_norm = 1e-8;
GN_options.TOL_t_k = 1e-7;
GN_options.max_iterations = 1000;
GN_options.t_k_start = 0.3;

[p_optim_end, optim_output] = ...
    GN_ass(F1_func, ...
           F2_func, ...
           optimization.p_optim_start(optimization.p_optim_estimable), ...
           optimization.lb(optimization.p_optim_estimable), ...
           optimization.ub(optimization.p_optim_estimable), ...
           GN_options);



