% forward integration with NURBS

L1 = 25.;
L2 = 0.;
L3 = 1.;

N3 = 50;

T_0 = 30;
T_end = 180;

heat_rate = 10.; % K/min

lambda_test_setup = [23*1, 35.6000, 0.9600];

k_sap_fit = [0.    45.5951   0.];

% c_p parametrization of sample
nrb_order = 4; % nrb_order = 4 equates to C^2

cntrl_pts = [10, 30, 60, 90, 120, 125, 130, 131, 135, 150, 160, 200; ...
             1, 1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52, 1.53, 1.58];
num_cntrl_pts = size(cntrl_pts,2);

% equidistant knots
knots = [zeros(1,nrb_order), ...
         (1:num_cntrl_pts-nrb_order)/(num_cntrl_pts-nrb_order+1), ...
         ones(1,nrb_order)];

c_p_sample = {'NURBS', [size(cntrl_pts,2), length(knots)]};

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup, ...
               'c_p_sample', c_p_sample};
p_sim = get_param_sim(common_args{:});

p_optim_start = [cntrl_pts(1,:), cntrl_pts(2,:), knots, k_sap_fit];

p_sim = update_c_p(p_sim, p_optim_start);

% T = 30:0.01:160;
% plot(T, p_sim.eval_c_p(T)); hold on
% plot(T, p_sim.eval_dc_p(T));



tic;
T_pcm = simulate_1d(p_sim);
toc


% Analytical reference solution
lambda_const = p_sim.lambda_test_setup(1);
rho_const = p_sim.rho_test_setup(1);
c_p_const = p_sim.c_p_test_setup(1);

heat_rate_s = heat_rate / 60;
dt = 0.05 / heat_rate_s; % fct evaluation every 0.05K
t = 0:dt:1/heat_rate_s*(T_end - T_0);
n = 100;
a = lambda_const / (c_p_const * rho_const);
T_ref = analytical_sol(L1,t,n,T_0, heat_rate_s, a); 



plot(T_ref,T_ref-T_pcm(:,p_sim.N1)); hold on


