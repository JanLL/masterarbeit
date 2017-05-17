function dT = ode_system1d(t, T, N1, N2, dx, heat_rate, lambda, J_lin_sparse)
% [dT] = ode_system1d(t, T, N, dx, heat_rate, lambda, J_lin_sparse)
% 
% Computes the right hand side of the 1D differential heat equation for
% density and specific heat capacity is temperature dependent and thermal
% conductivity is constant.
%
% INPUT:     t --> time
%            T --> temperature in degree Celsius.
%           N1 --> number of spatial discretization lattice points (part1).
%           N2 --> number of spatial discretization lattice points (PCM).
%           dx --> length [mm] of one spatial lattice point.
%    heat_rate --> rate [K/s] the temperature of the oven is increasing.
%       lambda --> thermal conductivity [mJ/mg*K]
% J_lin_sparse --> sparse matrix for the linear part pre-computed with
%                  build_linear_matrix(N) to avoid building up repeatedly
%                  at each function call.
%
% OUTPUT:   dT --> right hand side of the 1D differential heat equation
%                  \nabla \left[ \frac{\lambda}{\rho c_p} \nabla T \right]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de

% initial definitions
N = N1+N2;

c_p = ones(N, 1);
c_p(1:N1) = 0.41; % [mJ/(mg*K], Constantan, src: Wikipedia
c_p(N1+1:end) = c_p_formula(T(N1+1:end));

dc_p = zeros(N, 1);
dc_p(N1+1:end) = dc_p_formula(T(N1+1:end));


rho = ones(N, 1);
rho(1:N1) = 8.9; % [mg/mm^3], Constantan, src: Wikipedia
rho(N1+1:end) = rho_formula(T(N1+1:end));

drho = zeros(N, 1);
drho(N1+1:end) = drho_formula(T(N1+1:end));


lambda = ones(N, 1) * lambda;
lambda(1:N1) = 23; % [mW/(mm*K)]

%% Non-linear part
dT_non_lin = zeros(N,1);

dT_non_lin(1) = heat_rate;

% forward differences in gradient
dT_non_lin(N1+2:N-1) = ...
    (-lambda(N1+2:N-1) ./ (rho(N1+2:N-1) .* c_p(N1+2:N-1).^2) .* dc_p(N1+2:N-1) ...
     -lambda(N1+2:N-1) ./ (rho(N1+2:N-1).^2 .* c_p(N1+2:N-1)) .* drho(N1+2:N-1)) ...
    .* ((T(N1+3:N) - T(N1+2:N-1)).^2 ...
    ./ dx(N1+2:N-1).^2);

dT_non_lin(N) = 0;


%% Linear part

% linear part vector
dT_lin = lambda(1:N) ./ dx.^2 ./ (c_p(1:N) .* rho(1:N)) .* (J_lin_sparse * T(1:N));

%% put linear and non-linear part together

dT = dT_non_lin + dT_lin;


return