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
N = 500;
L = 5.; % [mm]
lambda = 0.96; % [mw/(mm * K)]
heat_rate = 0.3; % [K/min]
T_0 = 70.; % [min]
T_end = 200.; % [min]

% check for input arguments and update variables where necessary
if hasOption(varargin, 'N'), N = getOption(varargin, 'N'); end;
if hasOption(varargin, 'L'), L = getOption(varargin, 'L'); end;
if hasOption(varargin, 'lambda'), lambda = getOption(varargin, 'lambda'); end;
if hasOption(varargin, 'heat_rate'), heat_rate = getOption(varargin, 'heat_rate'); end;
if hasOption(varargin, 'T_0'), T_0 = getOption(varargin, 'T_0'); end;
if hasOption(varargin, 'T_end'), T_end = getOption(varargin, 'T_end'); end;


heat_rate = heat_rate / 60.;  % [K/min] -> [K/s]
dx = L/N;

% integration initial values
T0 = T_0 .* ones(N,1);

t0 = 0.;
tf = (T_end - T0(1)) / heat_rate;  % integrate up to T_oven = 200 degree Celsius
t = linspace(t0, tf, (tf-t0)*1)'; 

% pre-compute constant sparse matrix from linear part
J_lin_sparse = build_linear_matrix(N);


cols = ones(N, 3);
Jpattern = spdiags(cols, [0,1,-1], N, N);

opts = odeset('reltol', 1.0d-8, 'abstol', 1.0d-12, 'Jpattern', Jpattern);

tic;
sol = ode15s(@(t, y) ode_system1d(t, y, N, dx, heat_rate, lambda, J_lin_sparse), t, T0, opts);
toc


T = deval(sol, t)';


%figure(1);
% plot(t, T(:,1) - T(:,1), 'g', 'DisplayName', 'T_0'); hold on
% plot(t, T(:,2) - T(:,1), 'r', 'DisplayName', 'T_1');
% plot(t, T(:,3) - T(:,1), 'y', 'DisplayName', 'T_2');
%plot(T(:,1), T(:,end) - T(:,1), 'b', 'DisplayName', 'T_N'); hold on

%figure(2);
%plot(T(:,1), T(:,1)-T(:,2), 'DisplayName', 'T_0 - T_1'); hold on

% xlabel('Time t');
% ylabel('Temp T');
% legend(gca, 'show', 'Location', 'northwest')

if (nargout >= 1), varargout{1} = T; end
if (nargout >= 2), varargout{2} = sol; end

end

