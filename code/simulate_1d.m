function varargout = simulate_1d(eval_c_p, eval_dc_p, p_sim)
% [T, sol] = simulate_1d(varargin)
% 
% Solves the 1D heat equation for density rho and specific heat capacity
% temperature dependend. Spatial discretization is done via line method. 
%
% INPUT (needed):
%     eval_c_p --> fhandle to evaluate specific heat capacity.
%    eval_dc_p --> fhandle to evaluate derivative of specific heat capacity
%                  w.r.t. temperature.
%        p_sim --> struct with simulation parameters, available via the
%                  function get_default_sim_params()
%      
% OUTPUT:      T --> array of temperatures for all lattice points
% (variable) sol --> struct of solution of internal matlab ode solver.
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de

% export simulation parameters from struct to avoid writing p_sim.(...)
L1 = p_sim.L1;
L2 = p_sim.L2;
L3 = p_sim.L3;
N1 = p_sim.N1;
N2 = p_sim.N2;
N3 = p_sim.N3;
heat_rate = p_sim.heat_rate;
T_0 = p_sim.T_0;
T_end = p_sim.T_end;
c_p_test_setup = p_sim.c_p_test_setup;
rho_test_setup = p_sim.rho_test_setup;
lambda_test_setup = p_sim.lambda_test_setup;

heat_rate = heat_rate / 60.;  % [K/min] -> [K/s]
N = N1+N2+N3;
dx = ones(N, 1);
dx(1:N1)        = L1/N1;
dx(N1+1:N1+N2)  = L2/N2;
dx(N1+N2+1:end) = L3/N3;

% integration initial values
T0 = T_0 .* ones(N,1);
t0 = 0.;

tf = (T_end - T0(1)) / heat_rate;  % integrate up to T_oven = T_end degree Celsius
t = linspace(t0, tf, int32((T_end - T_0(1))*20. + 1))';
% function evaluation every 0.05 K of T_oven, independend of heat_rate

cols = ones(N, 3);
Jpattern = spdiags(cols, [0,1,-1], N, N);

opts = odeset('reltol', 1e-7, 'abstol', 1e-12, 'Jpattern', Jpattern);

% here changed either for PCM or Saphir!
ode_system1d_expl = @(t, y) ode_system1d_sap(t, y, N1, N2, N3, dx, heat_rate, ...
    c_p_test_setup, rho_test_setup, lambda_test_setup, eval_c_p, eval_dc_p);

%tic;
sol = ode15s(ode_system1d_expl, t, T0, opts);
%toc

T = deval(sol, t)';

if (nargout >= 1), varargout{1} = T; end
if (nargout >= 2), varargout{2} = sol; end

end

