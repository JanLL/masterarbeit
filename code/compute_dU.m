function [dU_output] = compute_dU(fit_data_input)
%
%
% INPUT:
%    fit_data_input --> either struct with saved fit data or string with
%                       path to .mat file with saved fit data.
%
% OUTPUT:
%                dU --> (n x 3) array where n number of T_ref evaluation
%                       points. 
%                       1st column: T_ref
%                       2nd column: dU from measurement
%                       3rd column: dU from optimization


if isstruct(fit_data_input)
    fit_data = fit_data_input;
elseif isstr(fit_data_input)
    fit_data = load(fit_data_input);
end

% abbreviations
dsc_data = fit_data.dsc_measurements.dsc_data;
index_T_dsc = fit_data.dsc_measurements.index_T_dsc;

p_sim = fit_data.simulation;
heat_rate = p_sim.heat_rate / 60; % [K/min] -> [K/s]
T_0 = p_sim.T_0;
T_end = p_sim.T_end;
lambda_const = p_sim.lambda_test_setup(1);
rho_const = p_sim.rho_test_setup(1);
c_p_const = p_sim.c_p_test_setup(1);
L1 = p_sim.L1;
N1 = p_sim.N1;

p_optim_end = fit_data.optimization.param_end;


% dU computation analog to compute_residuum (see there for details), 
% necessary for dU(T_ref) plot
dt = 0.05 / heat_rate; % fct evaluation every 0.05K
t = 0:dt:1/heat_rate*(T_end - T_0);
n = 100;
a = lambda_const / (c_p_const * rho_const);
T_ref = analytical_sol(L1,t,n,T_0, heat_rate, a);  

T_pcm = simulate_1d(p_sim.eval_c_p, p_sim.eval_dc_p, p_sim);

dT = T_ref(:) - T_pcm(:,N1);

k = p_sim.get_param_k(p_optim_end);
dU = polyval(k, dT);

index_T_p5 = find(T_ref(:) > p_sim.T_0 + 5, 1, 'first');

U_dsc = [dsc_data.data(index_T_dsc(1):index_T_dsc(2),1), dsc_data.data(index_T_dsc(1):index_T_dsc(2),3)];
dU_interp = interp1(T_ref(index_T_p5:end), dU(index_T_p5:end), ...
                    U_dsc(:,1), 'linear');

dU_output = [U_dsc(:,1), U_dsc(:,2), dU_interp];


end

