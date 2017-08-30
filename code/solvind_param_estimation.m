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
N1 = 200;
N3 = 50;  % error if N3=0
 
lambda_Const = 23.;  % [mW/(mm*K)]
rho_Const = 8.9;     % [mg/mm^3]
c_p_Const = 0.41;    % [mJ/(mg*K)]

lambda_pcm = 0.96;   % [mW/(mm*K)]
rho_pcm = 0.8;       % [mg/mm^3]

heat_rate = 10.;     % [K/min]
heat_rate_s = heat_rate / 60; % [K/min] -> [K/s]


% c_p parametrization with NURBS
% cntrl_pts_x = [0, 30, 60, 90, 120, 125, 130, 132., 135, 150, 160, 180];
% cntrl_pts_y = [1., 1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52, 1.53, 1.55];


% old best result for 10K/min, doesnt work in c++ step std::vector<double> C_temp = nurbs.eval_nurbs_curve(u);
cntrl_pts_x = [0,30,60,90,100:2:150,160,180,200];
cntrl_pts_y = [1.79858157644537,1.73123570378148,1.88939296848618,2.38582410198609,2.83324784260016,2.42651409950467,4.16171664906883,1.13860062268564,5.04919663138140,2.07976511011705,4.83068638838580,3.63392796870435,4.55952727527380,7.61255409096836,2.54742792701244,38.5762115714479,41.9573364402785,11.0500960538204,7.73154658124738,2.65644380692096,2.24473765106880,1.50806399242812,1.38297861923093,1.19283792707801,1.18634381049379,1.21363453226050,1.33993164075706,1.62786023313508,1.76526003655576,2.19534568189373,2.56165109798611,9.17894802559609,2];

num_cntrl_pts = length(cntrl_pts_x);


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
p_optim_start = [cntrl_pts_x, cntrl_pts_y];

% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_estimable(1:num_cntrl_pts) = false; % fix x-position of control pts

p_optim_fixed = p_optim_start(~p_optim_estimable);


%%%%%%%%%%% SolvIND initialization %%%%%%%%%%%%%%%%%%%%%%%
solvind('importDynamicModelLib', '/home/argo/SOLVIND_SUITE/Packages/SOLVIND/Debug/TEST/MODELS/libdynModelDesc_heat1D_pcm.so');

% Create model with grid specified by L1, L3, N1 and N3.
% Afterwards create integrator and link integrator with model.
model = solvind('createDynamicModel', 'heat1D_pcm', ...
                sprintf('0 %2.2f %2.2f %d %d %d -', L1, L3, N1, N3, num_cntrl_pts));
int = solvind('createIntegrator', 'daesol2_sparse_withCorrIters');
solvind('setModel', int, model);

model_dims = solvind('getDims', model);
num_params = model_dims.np;

% options
solvind('setTapeStorageMode', int, 'values');
solvind('setPrintLevel', int, 0);
pL = solvind('getPrintLevel', int);

solvind('setRelTol', int, 1e-3);  
% abs. Differenz (der integrierten FUnktion) von 1e-2 zu 1e-3 weniger als ein Promille
solvind('setMaxIntSteps', int, 2000);
solvind('setCorrectorAccuracyFactor', int, 1e-4);
solvind('setCorrectorAbsoluteAccuracy', int, 1e-10);

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
solvind_compute_residuum_expl(p_optim_start(p_optim_estimable));
solvind('reset');
return


% tape = solvind('getTape', int, 0);

% free memory
solvind('reset');

