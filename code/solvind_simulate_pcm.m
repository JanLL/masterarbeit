more off

% Set Simulation parameters
L1 = 15;  % [mm]
L3 = 0.5;  % [mm]
N1 = 750;
N3 = 50;

lambda_Const = 23.;  % [mW/(mm*K)]
rho_Const = 8.9;     % [mg/mm^3]
c_p_Const = 0.41;    % [mJ/(mg*K)]

lambda_pcm = 0.96;   % [mW/(mm*K)]
rho_pcm = 0.8;       % [mg/mm^3]

heat_rate = 10.;     % [K/min]

% c_p parametrization with NURBS
cntrl_pts_x = [0, 30, 60, 90, 120, 125, 130, 132., 135, 150, 160, 180];
cntrl_pts_y = [1., 1,  1.1, 1.15, 1.2, 5., 10, 1.5, 1.51, 1.52, 1.53, 1.55];
num_cntrl_pts = length(cntrl_pts_x);

% some pre-calculations
N = N1+N3;
a_Const = lambda_Const/(rho_Const*c_p_Const);


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
solvind('setMaxIntSteps', int, 100);
solvind('setCorrectorAccuracyFactor', int, 1e-3);
solvind('setCorrectorAbsoluteAccuracy', int, 1e-13);

solvind('setTimeHorizon', int, [0, 10.]);

ogrid = linspace(0, 10., 1001);
solvind('setContOutputConfig', int, ogrid);
solvind('storeAdjSensAtGrid', int);

T_0 = 30*ones(N,1);
p = [a_Const, lambda_pcm, rho_pcm, heat_rate, cntrl_pts_x, cntrl_pts_y]';

assert(num_params == length(p), ...
    'Number of paramters inconsistent in model and matlab code!');

initValues = [T_0; p];
solvind('setInitVals', int, initValues);

retval = solvind('evaluate', int);
% 
% if retval == 0
% 	sol = solvind('getSolution', int);
% 	contsol = solvind('getContOutput', int);
% 
% 	stats = solvind('getStats', int);
% 	timings = solvind('getTimings', int);
% end


% % compute first order forward sensitivities
% fwdSensDir = [zeros(1,N+num_params); eye(N+num_params)];
% solvind('setForwardTaylorCoefficients', int, N+num_params, 1, fwdSensDir);
% retval = solvind('forwardSensSweep', int);
% if retval == 0
% 	fwdSens = solvind('getFwdSens', int);
% end


% compute first order adjoints
% adjSensDir = eye(N);
% solvind('setForwardTaylorCoefficients', int, []);
% solvind('setAdjointTaylorCoefficients', int, N, 0, 1, adjSensDir);
% retval = solvind('backwardSensSweep', int);
% if retval == 0
% 	adjSens = solvind('getAdjSens', int);
% 	[gridT, gridAdj] = solvind('getAdjSensAtGrid', int);
% end





% fprintf('Deviation of fwd and adj 1st order derivatives: %e\n', ...
% 	norm(adjSens(2:end,:)' - fwdSens));

% tape = solvind('getTape', int, 0);


% free memory
solvind('reset');

