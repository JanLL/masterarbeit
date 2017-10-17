%%%%%%%%%% Get measurement data %%%%%%%%%%%%%%%%%%%
dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
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
L1 = 15;  % [mm]
L3 = 0.5;  % [mm]
N1 = 500;
N3 = 50;  % error if N3=0
 
lambda_Const = 23.;  % [mW/(mm*K)]
rho_Const = 8.9;     % [mg/mm^3]
c_p_Const = 0.41;    % [mJ/(mg*K)]

lambda_pcm = 0.96;   % [mW/(mm*K)]
rho_pcm = 0.85;      % [mg/mm^3]

heat_rate = 20.;     % [K/min]
T_0 = 10.;           % Start temperature oven [degC]
T_end = 200.;        % End temperature oven [degC]

sim_params_vec = [L1, L3, N1, N3, lambda_Const, rho_Const, c_p_Const, ...
                  lambda_pcm, m_pcm, ...
                  heat_rate, T_0, T_end];

heat_rate_s = heat_rate / 60; % [K/min] -> [K/s]

% c_p parametrization with Fraser-Suzuki-Peak
h  =  10.0;
r  =  100.0;
wr =  15.0;
sr =   0.1;
z  = 100.0;
b  =   2.0;
p_fraser_suzuki = [h, r, wr, sr, z, b];

p_atan_cp = [125., 10, 0.01, 10., 0.003, 2];

p_gauss_lin_comb = [10, 0.1, 130, ...
                    1,  1,   128, ...
                    0.5,  1,   126, ...
                    0.1,  1,   124, ...
                    0.05, 1,   122, ...
                    0., 2.];
                
fit_data = load('/home/argo/masterarbeit/fits_data/2017-10-15_23:26:52_407_0,3Kmin_L1=50_L3=0,1/fit_data.mat');
p_gauss_lin_comb = fit_data.optimization.p_optim_end;
% p_gauss_lin_comb = [38.2793,    1.2498,  123.0044, ...
%                     1.1022,   13.2169,  104.8824, ...
%                     0.,    3.7355,  124.3034, ...
%                     3.8868, 7.1009,  118.2194, ...
%                     0.5846,   30.5770,   76.8653, ...
%                     0., 1., 0., ...
%                     0., 1., 0., ...
%                     0., 1., 0., ...
%                     0., 1., 0., ...
%                     0., 1., 0., ...
%                     1., 1.6483];


%%%%%%%%%% some pre-calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = N1+N3;
dx_Const = L1/N1;
dx_pcm = L3/N3;
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


% Set optimization variables
c_p_param_type = 'gauss_linear_comb';
p_optim_start = p_gauss_lin_comb;

% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_fixed = p_optim_start(~p_optim_estimable);

num_free_optim_params = sum(p_optim_estimable);



heat1D_pcm('reset');
heat1D_pcm('init', sim_params_vec, meas_data, c_p_param_type);


figure(1); % q_pcm_in plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();
figure(3); % dqdp plot
ax3 = gca();

compute_q_dqdp_mex_expl = @(p_optim) compute_q_dqdp_mex(...
    p_optim, p_optim_estimable, p_optim_fixed, c_p_param_type, T_ref_dsc, q_dsc, ax1, ax2, ax3);


%%%%%%%%%%%%%% INITIAL VALUE TEST %%%%%%%%%%%%%%%%%%%%
[res, dqdp] = compute_q_dqdp_mex_expl(p_optim_start);

q_sim_500 = res + q_dsc;

return

% close all;
% relErr_dqdp = abs(1 - dqdp.fwd ./ dqdp.adj);
% relErr_dqdp(isnan(relErr_dqdp)) = 0.;
% figure()
% image(relErr_dqdp, 'CDataMapping', 'scaled')
% colorbar

close all
figure()
sens_diff = load('/home/argo/masterarbeit/diff_sens.txt');
max(max(abs(sens_diff)))
%imagesc(sens_diff, [0, 1e-12]);
image(sens_diff, 'CDataMapping', 'scaled')
colorbar


sens = load('/home/argo/masterarbeit/sens.txt');
fwdSens_N1 = sens(:,1:4:end);
adjSens_N1 = sens(:,2:4:end);
fwdSens_N1p1 = sens(:,3:4:end);
adjSens_N1p1 = sens(:,4:4:end);

figure()
image(fwdSens_N1, 'CDataMapping', 'scaled'); colorbar
figure()
image(adjSens_N1, 'CDataMapping', 'scaled'); colorbar
figure()
image(fwdSens_N1p1, 'CDataMapping', 'scaled'); colorbar
figure()
image(adjSens_N1p1, 'CDataMapping', 'scaled'); colorbar

fwdSens_N1(772,2) / adjSens_N1(772,2) - 1


return;

% Jacobian via finite differences
% [x,~,~,~,~,~,dqdp_finite_diff] = lsqnonlin(...
%     compute_q_dqdp_mex_expl, p_optim_start(p_optim_estimable), [], [], optimset('MaxIter',0));
% figure()
% image(dqdp_finite_diff, 'CDataMapping', 'scaled');
% colorbar
% return;


%%%%%%%%%%%%%%%%%%%%%%% SOLVE OPTIMIZATION PROBLEM %%%%%%%%%%%%%%%%%%%%%%%
lb = ones(1,num_free_optim_params) * -1;
ub = ones(1,num_free_optim_params) * 500.;


%ub(4) = 1.;  % sr < 1
optim_con = {lb, ub};

opt_options = optimoptions('lsqnonlin', ...
                           'Display', 'iter-detailed', ...
                           'OutputFcn', @disp_aux, ...
                           'SpecifyObjectiveGradient', true);
[p_optim,~,~,~,optim_output,~,jac_output] = lsqnonlin(...
    compute_q_dqdp_mex_expl, p_optim_start(p_optim_estimable), optim_con{:}, opt_options);




