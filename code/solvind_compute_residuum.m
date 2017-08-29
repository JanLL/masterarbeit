function [residuum, dq_dp] = solvind_compute_residuum(...
    p_optim_free, p_optim_estimable, p_optim_fixed, p_sim, int, q_dsc, ...
    m_pcm, ax1, ax2)

% TODO: Description!

% Abbreviations
L1 = p_sim.L1;
L3 = p_sim.L3;
N1 = p_sim.N1;
N3 = p_sim.N3;
a_Const = p_sim.a_Const;
lambda_pcm = p_sim.lambda_pcm;
rho_pcm = p_sim.rho_pcm;
m_pcm = p_sim.m_pcm;
heat_rate = p_sim.heat_rate;
T_0 = p_sim.T_0;

N = N1+N3;
dx_Const = L1/N1;
dx_pcm = L3/N3;
num_meas = length(q_dsc(:,1));

% build vector of all (free and fixed) optimization parameters
p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim_free;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

num_cntrl_pts = length(p_optim_all) / 2;
cntrl_pts_x = p_optim_all(1:num_cntrl_pts);
cntrl_pts_y = p_optim_all(num_cntrl_pts+1:end);

T_0_vec = T_0*ones(N,1); 
p = [a_Const, lambda_pcm, rho_pcm, heat_rate, cntrl_pts_x, cntrl_pts_y]';
num_params_all = length(p);

initValues = [T_0_vec; p];
solvind('setInitVals', int, initValues);

fprintf('Performing forward integration.\n');
retval = solvind('evaluate', int);

if retval == 0
	%sol = solvind('getSolution', int);
	T_sol = solvind('getContOutput', int);

	%stats = solvind('getStats', int);
	%timings = solvind('getTimings', int);
end

constant_factor_heat_flux = (lambda_pcm*m_pcm)/(rho_pcm*dx_pcm*N3);

% Compute Residuum
q_sim = -constant_factor_heat_flux * (T_sol(N1+2,:) - T_sol(N1+1,:));
residuum = q_dsc(:,2)' - q_sim;


% Compute first order forward sensitivities
fprintf('Computing first order forward sensitivities.\n');
fwdSensDir = [zeros(1,N+num_params_all); eye(N+num_params_all)];
solvind('setForwardTaylorCoefficients', int, N+num_params_all, 1, fwdSensDir);
retval = solvind('forwardSensSweep', int);
if retval == 0
	fwdSens = solvind('getFwdSens', int);
end

dT_dp = fwdSens(:,N1+N3+1:end); 

% Build Jacobian of resiidum w.r.t. free parameters
row_index = [1:num_meas, 1:num_meas];
col_index = [(N1+1)*ones(1, num_meas), (N1+2)*ones(1, num_meas)];
values = [ones(1, num_meas), -1.*ones(1, num_meas)] * constant_factor_heat_flux;
dq_dT_sparse = sparse(row_index, col_index, values, num_meas, N1+N3);

dq_dp = dq_dT_sparse * dT_dp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cla(ax1); % q_pcm_in plot
hold(ax1, 'on')
plot(ax1, q_dsc(:,1), q_sim, 'DisplayName', 'Simulation');
plot(ax1, q_dsc(:,1), q_dsc(:,2), 'DisplayName', 'Measurement');
plot(ax1, q_dsc(:,1), residuum, 'DisplayName', 'Residuum');
legend(ax1, 'show', 'location', 'northoutside');
xlabel(ax1, 'T_{ref} [degC]');
ylabel(ax1, 'q_{pcm}^{in} [mW]');
drawnow;




end