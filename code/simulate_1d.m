function varargout = simulate_1d(eval_c_p, eval_dc_p, varargin)
% [T, sol] = simulate_1d(varargin)
% 
% Solves the 1D heat equation for density rho and specific heat capacity
% temperature dependend. Spatial discretization is done via line method. 
%
% INPUT (needed):
%     eval_c_p --> fhandle to evaluate specific heat capacity.
%    eval_dc_p --> fhandle to evaluate derivative of specific heat capacity
%                  w.r.t. temperature.%
% INPUT (variable):          
%        N1 --> number of spatial discretization lattice points (Constantan)
%        N2 --> number of spatial discretization lattice points (Crucible)
%        N3 --> number of spatial discretization lattice points (PCM)
%        L1 --> length [mm] of Constantan part.
%        L2 --> length [mm] of crucible part.
%        L3 --> length [mm] of PCM part.
%       T_0 --> initial temperature [degree celsius] everywhere.
%     T_end --> integrate until T_oven as reached T_end [degree celsius].
% heat_rate --> rate [K/min] with which oven temperature increases.
%      
% OUTPUT:      T --> array of temperatures for all lattice points
% (variable) sol --> struct of solution of internal matlab ode solver.
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de


% check for input arguments and update variables where necessary
if hasOption(varargin, 'L1'), L1 = getOption(varargin, 'L1'); 
else L1 = 25.; end
if hasOption(varargin, 'L2'), L2 = getOption(varargin, 'L2'); 
else L2 = 0.; end
if hasOption(varargin, 'L3'), L3 = getOption(varargin, 'L3'); 
else L3 = 1.; end

if hasOption(varargin, 'N3'), N3 = getOption(varargin, 'N3');
else N3 = 50; end

% compute N1, N2 s.t. dx is equal everywhere as default
if hasOption(varargin, 'N1'), N1 = getOption(varargin, 'N1'); 
else N2 = L2 / L3 * N3; end
if hasOption(varargin, 'N2'), N2 = getOption(varargin, 'N2'); 
else N1 = L1 / L3 * N3; end

if hasOption(varargin, 'heat_rate'), heat_rate = getOption(varargin, 'heat_rate'); 
else heat_rate = 10; end
if hasOption(varargin, 'T_0'), T_0 = getOption(varargin, 'T_0'); 
else T_0 = 10.; end
if hasOption(varargin, 'T_end'), T_end = getOption(varargin, 'T_end'); 
else T_end = 300.; end


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
t = linspace(t0, tf, int32((T_end - T_0(1))*20.))';
% function evaluation every 0.05 K, independend of heat_rate

cols = ones(N, 3);
Jpattern = spdiags(cols, [0,1,-1], N, N);

opts = odeset('reltol', 1e-10, 'abstol', 1e-12, 'Jpattern', Jpattern);

ode_system1d_expl = @(t, y) ode_system1d(t, y, N1, N2, N3, dx, heat_rate, ...
    eval_c_p, eval_dc_p);

tic;
sol = ode15s(ode_system1d_expl, t, T0, opts);
toc

T = deval(sol, t)';

if (nargout >= 1), varargout{1} = T; end
if (nargout >= 2), varargout{2} = sol; end

end

