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
T_0 = 10.;
T_end = 200.;

L1 = 15;  % [mm]
L3 = 0.5;  % [mm]
N1 = 300;
N3 = 50;  % error if N3=0
 
lambda_Const = 23.;  % [mW/(mm*K)]
rho_Const = 8.9;     % [mg/mm^3]
c_p_Const = 0.41;    % [mJ/(mg*K)]

lambda_pcm = 0.96;   % [mW/(mm*K)]
rho_pcm = 0.8;       % [mg/mm^3]

heat_rate = 10.;     % [K/min]
heat_rate_s = heat_rate / 60; % [K/min] -> [K/s]

% c_p parametrization with Fraser-Suzuki-Peak
h  =  50.0;
r  =  35.0;
wr =  15.0;
sr =   0.2;
z  = 125.0;
b  =   2.0;

p_fraser_suzuki = [h, r, wr, sr, z, b];

% c_p parametrization with old atan function
p_atan_cp = [125., 10, 0.01, 10., 0.003, 2];



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

heat1D_pcm('init', meas_data, 'fraser-suzuki');



figure(1); % q_pcm_in plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();

% Set optimization variables
p_optim_start = p_fraser_suzuki;

% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_fixed = p_optim_start(~p_optim_estimable);

num_free_optim_params = sum(p_optim_estimable);


compute_q_dqdp_mex_expl = @(p_optim) compute_q_dqdp_mex(...
    p_optim, p_optim_estimable, p_optim_fixed, T_ref_dsc, q_dsc, ax1, ax2);


%%%% INITIAL VALUE TEST %%%%%%%%%%%%
compute_q_dqdp_mex_expl(p_optim_start);
return;





