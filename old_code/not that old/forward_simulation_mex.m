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
L1 = 40;  % [mm]
L3 = 0.1;  % [mm]

n_pcm = 0.2;
N = 350;
N3 = N*n_pcm;
N1 = N - N3;

% N3 = 50;  % error if N3=0
% N1 = 2500;
% N = N1 + N3;

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
n_tr = 0.03;
n_m = 0.01;
t = 1.;

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

sum(spatial_gridsize(1:N1-1))
sum(spatial_gridsize(1:N1+N3-1))
figure(10)
plot(spatial_gridsize, 'x')
return


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
                
fit_data = load(['/home/argo/masterarbeit/fits_data/', ...
                '2017-10-15_22:45:00_407_L1=50_L3=0.1_N1=200_N3=50/', ...
                '2017-10-15_23:26:52_407_0,3Kmin_L1=50_L3=0,1/fit_data.mat']);

h_tilde = 1.;            
p_gauss_lin_comb = [fit_data.optimization.p_optim_end, h_tilde];
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
% N = N1+N3;
% dx_Const = L1/N1;
% dx_pcm = L3/N3;
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
p_optim_free = p_optim_start(p_optim_estimable);

num_free_optim_params = sum(p_optim_estimable);



heat1D_pcm('reset');
heat1D_pcm('init', sim_params_vec, spatial_gridsize, meas_data, c_p_param_type);


figure(1); % q_pcm_in plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();
figure(3); % dqdp plot
ax3 = gca();



p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim_free;
p_optim_all(~p_optim_estimable) = p_optim_fixed;


% c_p plot
T_plot = 30:0.01:160;
switch c_p_param_type
    case 'old_atan_formula'
        c_p_plot = c_p_formula(T_plot, p_optim_all(1:6));
    case 'fraser_suzuki'
        c_p_plot = c_p_fs(T_plot, p_optim_all);
    case 'gauss_linear_comb'
        c_p_plot = c_p_gauss_linear_comb(T_plot, p_optim_all);
end

cla(ax2); 
hold(ax2, 'on')
plot(ax2, T_plot, c_p_plot, 'DisplayName', 'c_p Simulation');
legend(ax2, 'show', 'location', 'northwest');
xlabel(ax2, 'T [degC]');
ylabel(ax2, 'c_p [mJ/(mg*K]');
drawnow;

%[residuum, Jac] = heat1D_pcm('optimization', p_optim_all); 
[residuum, Jac, T] = heat1D_pcm('optimization', p_optim_all); 
Jac = Jac(:,p_optim_estimable);

T_350 = T(:,N1);

%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%
% heat flux plot
q_sim = residuum + q_dsc;
cla(ax1); % q_pcm_in plot
hold(ax1, 'on')
plot(ax1, T_ref_dsc, q_sim, 'DisplayName', 'Simulation');
plot(ax1, T_ref_dsc, q_dsc, 'DisplayName', 'Measurement');
plot(ax1, T_ref_dsc, residuum, 'DisplayName', 'Residuum');
legend(ax1, 'show', 'location', 'northwest');
xlabel(ax1, 'T_{ref} [degC]');
ylabel(ax1, 'q_{pcm}^{in} [mW]');
drawnow;


% dqdp plot
cla(ax3); 
hold(ax3, 'on')
image(ax3, Jac, 'CDataMapping', 'scaled');
colorbar(ax3);
title(ax3, 'dqdp')
xlabel(ax3, 'c_p parameters');
ylabel(ax3, 'T_{ref}');
drawnow;



%% Compute relative error of different grid sizes
close all;
figure()
plot(abs(1 - T_350 ./ T_2500), 'x')
set(gca,'FontSize',10);
xlabel('measurement points');
ylabel('Absolute relative error of T^{N_1}');
title(sprintf('N1=350,2500, N3=%d, n_{tr}=%1.2f, n_m=%1.2f, t=%1.3f',N3,n_tr,n_m,t));








