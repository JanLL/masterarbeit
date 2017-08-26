more off

% Set Simulation parameters
L1 = 15; % [mm]
L3 = 0.5; % [mm]
N1 = 201;
N3 = 50;



solvind('importDynamicModelLib', '/home/argo/SOLVIND_SUITE/Packages/SOLVIND/Debug/TEST/MODELS/libdynModelDesc_heat1D_pcm.so');

% grid sizes for level 0
model = solvind('createDynamicModel', 'heat1D_pcm', sprintf('0 %d %d %d %d -', L1, L3, N1, N3));

int = solvind('createIntegrator', 'daesol2_sparse_withCorrIters');


dims = solvind('getDims', int);
num_params = dims.np;


solvind('setModel', int, model);
solvind('setTapeStorageMode', int, 'values');
solvind('setPrintLevel', int, 0);
pL = solvind('getPrintLevel', int);

solvind('setRelTol', int, 1e-3);
solvind('setMaxIntSteps', int, 100);
solvind('setCorrectorAccuracyFactor', int, 1e-3);
solvind('setCorrectorAbsoluteAccuracy', int, 1e-13);

solvind('setTimeHorizon', int, [0 0.1]);

ogrid = linspace(0, 0.1, 21);
solvind('setContOutputConfig', int, ogrid);
solvind('storeAdjSensAtGrid', int);

y0 = [0.99*ones(N,1); ones(2,1)];
solvind('setInitVals', int, y0);

retval = solvind('evaluate', int);

if retval == 0
	sol = solvind('getSolution', int);

	contsol = solvind('getContOutput', int);
	%[X,Y] = meshgrid(ogrid, linspace(0, 1, 11));
	%surf(X, Y, contsol);
	%drawnow

	stats = solvind('getStats', int);
	timings = solvind('getTimings', int);
end



% compute first order forward sensitivities
fwdSensDir = [zeros(1,N+2); eye(N+2)];
solvind('setForwardTaylorCoefficients', int, N+2, 1, fwdSensDir);
retval = solvind('forwardSensSweep', int);
if retval == 0
	fwdSens = solvind('getFwdSens', int);
end


% compute first order adjoints
adjSensDir = eye(N);
solvind('setForwardTaylorCoefficients', int, []);
solvind('setAdjointTaylorCoefficients', int, N, 0, 1, adjSensDir);
retval = solvind('backwardSensSweep', int);
if retval == 0
	adjSens = solvind('getAdjSens', int);
	[gridT, gridAdj] = solvind('getAdjSensAtGrid', int);
end





fprintf('Deviation of fwd and adj 1st order derivatives: %e\n', ...
	norm(adjSens(2:end,:)' - fwdSens));

tape = solvind('getTape', int, 0);


% free memory
solvind('reset');

