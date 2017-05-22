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

% Default values
% Constantan
N1 = 50;
L1 = 50.; % [mm]

% crucible
N2 = 100.;
L2 = 2.; % [mm]

% PCM
N3 = 750;
L3 = 5.; % [mm]

lambda = 0.96; % [mW/(mm * K)]
heat_rate = 10.; % [K/min]
T_0 = 10.; % [degree Celsius]
T_end = 300.; % [degree Celsius]

% check for input arguments and update variables where necessary
if hasOption(varargin, 'N1'), N1 = getOption(varargin, 'N1'); end;
if hasOption(varargin, 'L1'), L1 = getOption(varargin, 'L1'); end;
if hasOption(varargin, 'N2'), N2 = getOption(varargin, 'N2'); end;
if hasOption(varargin, 'L2'), L2 = getOption(varargin, 'L2'); end;
if hasOption(varargin, 'N3'), N3 = getOption(varargin, 'N3'); end;
if hasOption(varargin, 'L3'), L3 = getOption(varargin, 'L3'); end;
if hasOption(varargin, 'lambda'), lambda = getOption(varargin, 'lambda'); end;
if hasOption(varargin, 'heat_rate'), heat_rate = getOption(varargin, 'heat_rate'); end;
if hasOption(varargin, 'T_0'), T_0 = getOption(varargin, 'T_0'); end;
if hasOption(varargin, 'T_end'), T_end = getOption(varargin, 'T_end'); end;


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
t = linspace(t0, tf, (tf-t0)*1)'; 

% pre-compute constant sparse matrix from linear part
J_lin_sparse = build_linear_matrix(N);


cols = ones(N, 3);
Jpattern = spdiags(cols, [0,1,-1], N, N);

opts = odeset('reltol', 1e-8, 'abstol', 1e-8, 'Jpattern', Jpattern);

tic;
sol = ode15s(@(t, y) ode_system1d(t, y, N1, N2, N3, dx, heat_rate, lambda, J_lin_sparse), t, T0, opts);
toc

T = deval(sol, t)';

%plot(t, T(:,N1));

if (nargout >= 1), varargout{1} = T; end
if (nargout >= 2), varargout{2} = sol; end

end

