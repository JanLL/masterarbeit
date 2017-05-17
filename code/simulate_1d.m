function varargout = simulate_1d(varargin)
% [T, sol] = simulate_1d(varargin)
% 
% Solves the 1D heat equation for density rho and specific heat capacity
% temperature dependend. Spatial discretization is done via line method. 
%
% INPUT:          N --> number of spatial discretization lattice points.
% (variable)      L --> length [mm] of line where the heat diffuses.
%            lambda --> thermal conductivity [mw/(mm * K)]
%            heat_rate --> rate [K/min] with which oven temperature increases
%
% OUTPUT:      T --> array of temperatures for all lattice points
% (variable) sol --> struct of solution of internal matlab ode solver.
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de


startup;  % load paths to external .m functions

% default values
N1 = 50;
L1 = 50.; % [mm]
N2 = 750;
L2 = 5.; % [mm]
lambda = 0.96; % [mw/(mm * K)]
heat_rate = 0.3; % [K/min]
T_0 = 10.; % [min]
T_end = 200.; % [min]

% check for input arguments and update variables where necessary
if hasOption(varargin, 'N1'), N1 = getOption(varargin, 'N1'); end;
if hasOption(varargin, 'L1'), L1 = getOption(varargin, 'L1'); end;
if hasOption(varargin, 'N2'), N2 = getOption(varargin, 'N2'); end;
if hasOption(varargin, 'L2'), L2 = getOption(varargin, 'L2'); end;
if hasOption(varargin, 'lambda'), lambda = getOption(varargin, 'lambda'); end;
if hasOption(varargin, 'heat_rate'), heat_rate = getOption(varargin, 'heat_rate'); end;
if hasOption(varargin, 'T_0'), T_0 = getOption(varargin, 'T_0'); end;
if hasOption(varargin, 'T_end'), T_end = getOption(varargin, 'T_end'); end;


heat_rate = heat_rate / 60.;  % [K/min] -> [K/s]
N = N1+N2;
dx = ones(N, 1);
dx(1:N1) = L1/N1;
dx(N1+1:end) = L2/N2;

% integration initial values
T0 = T_0 .* ones(N,1);

t0 = 0.;
tf = (T_end - T0(1)) / heat_rate;  % integrate up to T_oven = 200 degree Celsius
t = linspace(t0, tf, (tf-t0)*1)'; 

% pre-compute constant sparse matrix from linear part
J_lin_sparse = build_linear_matrix(N);


cols = ones(N, 3);
Jpattern = spdiags(cols, [0,1,-1], N, N);

opts = odeset('reltol', 1.0d-6, 'abstol', 1.0d-6, 'Jpattern', Jpattern);

tic;
sol = ode15s(@(t, y) ode_system1d(t, y, N1, N2, dx, heat_rate, lambda, J_lin_sparse), t, T0, opts);
toc


T = deval(sol, t)';

if (nargout >= 1), varargout{1} = T; end
if (nargout >= 2), varargout{2} = sol; end

end

