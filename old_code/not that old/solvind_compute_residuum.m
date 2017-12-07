function [residuum, dq_dp] = solvind_compute_residuum(...
    p_optim_free, p_optim_estimable, p_optim_fixed, p_sim, int, q_dsc, ax1, ax2)

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

num_optim_params = length(p_optim_estimable);

h = p_optim_all(1);
r = p_optim_all(2);
wr = p_optim_all(3);
sr = p_optim_all(4);
z = p_optim_all(5);
b = p_optim_all(6);


T_0_vec = T_0*ones(N,1); 
p = [a_Const, lambda_pcm, rho_pcm, heat_rate, h, r, wr,  sr, z, b]';

initValues = [T_0_vec; p];
solvind('setInitVals', int, initValues);

fprintf('Performing forward integration.\n');
retval = solvind('evaluate', int);

if retval == 0
    fprintf('Integration successful!\n');
	%sol = solvind('getSolution', int);
	T_sol = solvind('getContOutput', int);
    
	%stats = solvind('getStats', int);
	%timings = solvind('getTimings', int);
end

constant_factor_heat_flux = (lambda_pcm*m_pcm)/(rho_pcm*dx_pcm^2*N3);

% Compute Residuum
q_sim = -constant_factor_heat_flux * (T_sol(N1+2,:) - T_sol(N1+1,:));
residuum = q_dsc(:,2)' - q_sim;

if nargout > 1
    % Compute first order forward sensitivities
    fprintf('Computing first order forward sensitivities.\n');
    %fwdSensDir = [zeros(1,N+num_params_all); eye(N+num_params_all)];
    fwdSensDir = [zeros(1,num_optim_params); zeros(N+4, num_optim_params); diag(p_optim_estimable)];
    fwdSensDir(:, ~any(fwdSensDir,1)) = [];  % remove zero-columns

    solvind('setForwardTaylorCoefficients', int, size(fwdSensDir, 2), 1, fwdSensDir);
    tic;
    retval = solvind('forwardSensSweep', int);
    toc
    if retval == 0
        fwdSens = solvind('getFwdSens', int);
    end

    dT_dp = fwdSens(:,:); 

    % Build Jacobian of resiidum w.r.t. free parameters
    row_index = [1:num_meas, 1:num_meas];
    col_index = [(N1+1)*ones(1, num_meas), (N1+2)*ones(1, num_meas)];
    values = [ones(1, num_meas), -1.*ones(1, num_meas)] * constant_factor_heat_flux;
    dq_dT_sparse = sparse(row_index, col_index, values, num_meas, N1+N3);
    
    dq_dp = - dq_dT_sparse * dT_dp;  % minus because: residuum = q_dsc(:,2)' - q_sim;
        
%     figure()
%     image(dT_dp, 'CDataMapping', 'scaled')
%     colorbar
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heat flux plot
cla(ax1); % q_pcm_in plot
hold(ax1, 'on')
plot(ax1, q_dsc(:,1), q_sim, 'DisplayName', 'Simulation');
plot(ax1, q_dsc(:,1), q_dsc(:,2), 'DisplayName', 'Measurement');
plot(ax1, q_dsc(:,1), residuum, 'DisplayName', 'Residuum');
legend(ax1, 'show', 'location', 'northoutside');
xlabel(ax1, 'T_{ref} [degC]');
ylabel(ax1, 'q_{pcm}^{in} [mW]');
drawnow;

% c_p plot
T_plot = 30:0.01:160;
c_p_plot = frasersuzuki(T_plot, [h, r, z, wr, sr]) + b;
%c_p_plot = c_p_formula(T_plot, p_optim_all(1:6));

cla(ax2); 
hold(ax2, 'on')
plot(ax2, T_plot, c_p_plot, 'DisplayName', 'c_p Simulation');
legend(ax1, 'show', 'location', 'northoutside');
xlabel(ax2, 'T [degC]');
ylabel(ax2, 'c_p [mJ/(mg*K]');
drawnow;

end