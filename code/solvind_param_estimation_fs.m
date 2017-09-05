more off

%%%%%%%%%% Get measurement data %%%%%%%%%%%%%%%%%%%
dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
dsc = DSC204_readFile(dsc_filename);

m_pcm = dsc.mass;

% TODO: sinnvolles Intervall automatisch waehlen ... wobei das hier fuer
% alle Messungen bisher ganz gut war
index_T_dsc = [find(dsc.data(:,1) > 29, 1, 'first'), ...
               find(dsc.data(:,1) < 157.9, 1, 'last')];
q_dsc = [dsc.data(index_T_dsc(1):index_T_dsc(2),1), ...
         dsc.data(index_T_dsc(1):index_T_dsc(2),3) ...
         ./ dsc.data(index_T_dsc(1):index_T_dsc(2),4)];
q_dsc(:,2) = q_dsc(:,2) * m_pcm;

% q_dsc(:,1): measurement times
% q_dsc(:,2): heat flux values  

num_meas = length(q_dsc(:,1));


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

meas_times = zeros(1,num_meas);
for i=1:length(q_dsc(:,1))

    F = @(t) analytical_sol(L1,t,n,T_0, heat_rate_s, a) - q_dsc(i,1);

    t_guess = (q_dsc(i,1) - T_0)/heat_rate_s;
    
    fsolve_options = optimoptions('fsolve','Display','none');
    meas_times(i) = fsolve(F, t_guess, fsolve_options);
end


% Create simulation parameter struct
p_sim = struct();

p_sim.L1 = L1;
p_sim.L3 = L3;
p_sim.N1 = N1;
p_sim.N3 = N3;
p_sim.a_Const = a_Const;
p_sim.lambda_pcm = lambda_pcm;
p_sim.rho_pcm = rho_pcm;
p_sim.m_pcm = m_pcm;
p_sim.heat_rate = heat_rate;
p_sim.T_0 = T_0;

% Set optimization variables
p_optim_start = p_fraser_suzuki;

% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_fixed = p_optim_start(~p_optim_estimable);

num_free_optim_params = sum(p_optim_estimable);


%%%%%%%%%%% SolvIND initialization %%%%%%%%%%%%%%%%%%%%%%%
solvind('importDynamicModelLib', '/home/argo/SOLVIND_SUITE/Packages/SOLVIND/Debug/TEST/MODELS/libdynModelDesc_heat1D_pcm.so');

% Create model with grid specified by L1, L3, N1 and N3.
% Afterwards create integrator and link integrator with model.
model = solvind('createDynamicModel', 'heat1D_pcm', ...
                sprintf('0 %2.2f %2.2f %d %d -', L1, L3, N1, N3));
int = solvind('createIntegrator', 'daesol2_sparse_withCorrIters');
solvind('setModel', int, model);

model_dims = solvind('getDims', model);
num_params = model_dims.np;

% options
solvind('setTapeStorageMode', int, 'values');
solvind('setPrintLevel', int, 0);
pL = solvind('getPrintLevel', int);


solvind('setMaxBDFOrder', int, 4);
solvind('setRelTol', int, 1e-6);
solvind('setMaxIntSteps', int, 2000);
solvind('setCorrectorAccuracyFactor', int, 1e-5);
solvind('setCorrectorAbsoluteAccuracy', int, 1e-13);

t_0 = 0.;
t_end = (T_end - T_0) / heat_rate_s;  % vllt an meas_times(end) koppeln...

solvind('setTimeHorizon', int, [t_0, t_end]);
solvind('setContOutputConfig', int, meas_times);  
% just evaluate at times where also measurements were done to cpmpute residuum later

%solvind('storeAdjSensAtGrid', int); % benutzen erstmal nur fwdSens


figure(1); % q_pcm_in plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();

solvind_compute_residuum_expl = @(p_optim) ... 
    solvind_compute_residuum(p_optim, p_optim_estimable, ...
        p_optim_fixed, p_sim, int, q_dsc, ax1, ax2);

% TEST INITIAL VALUE
[res, Jac] = solvind_compute_residuum_expl(p_optim_start(p_optim_estimable));
solvind('reset');
return


%%%%%%%%%%%%%%%%%%%%%%% SOLVE OPTIMIZATION PROBLEM %%%%%%%%%%%%%%%%%%%%%%%
lb = zeros(1,num_free_optim_params);
ub = ones(1,num_free_optim_params)*300.;
%ub(4) = 1.;  % sr < 1
optim_con = {lb, ub};

opt_options = optimoptions('lsqnonlin', ...
                           'Display', 'iter-detailed', ...
                           'OutputFcn', @disp_aux, ...
                           'SpecifyObjectiveGradient', true);
[p_optim,~,~,~,optim_output,~,jac_output] = lsqnonlin(...
    solvind_compute_residuum_expl, p_optim_start(p_optim_estimable), optim_con{:}, opt_options);


% tape = solvind('getTape', int, 0);

% free memory
solvind('reset');

