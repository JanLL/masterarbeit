function dT = ode_system1d(t, T, N1, N2, N3, dx, heat_rate, ...
                           eval_c_p, eval_dc_p, J_lin_sparse)
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
%   c_p_params --> Parameter of function for specific heat capacity
% J_lin_sparse --> sparse matrix for the linear part pre-computed with
%                  build_linear_matrix(N) to avoid building up repeatedly
%                  at each function call.
%
% OUTPUT:   dT --> right hand side of the 1D differential heat equation
%                  \nabla \left[ \frac{\lambda}{\rho c_p} \nabla T \right]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de

% initial definitions
N = N1+N2+N3;

c_p = ones(N, 1);
c_p(1:N1) = 0.41; % [mJ/(mg*K], Constantan, src: Wikipedia
c_p(N1+1:N1+N2) = 0.99; % [mJ/(mg*K], Al2O3, src: www.pgo-online.com
c_p(N1+N2+1:end) = eval_c_p(T(N1+N2+1:end));

dc_p = zeros(N, 1);
dc_p(N1+N2+1:end) = eval_dc_p(T(N1+N2+1:end));


rho = ones(N, 1);
rho(1:N1) = 8.9; % [mg/mm^3], Constantan, src: Wikipedia
%rho(1:N1) = 40.; % Um zu zeigen, dass eine anpassung des querschnitts ueber
                 % die Dichte einen unterschied macht...

rho(N1+1:N1+N2) = 3.75; % [mg/mm^3], Al2O3, src: www.pgo-online.com
rho(N1+N2+1:end) = rho_formula(T(N1+N2+1:end));

drho = zeros(N, 1);
drho(N1+N2+1:end) = drho_formula(T(N1+N2+1:end));


lambda = ones(N, 1) * 0.96; % [mJ/mg*K], PCM, src: Robert: PCM_lambda.m
lambda(1:N1) = 23. * 4.; % [mW/(mm*K)], Constantan, src: Wikipedia
lambda(N1+1:N1+N2) = 35.6; % [mW/(mm*K)], Al2O3, src: Wikipedia

%% Non-linear part
dT_non_lin = zeros(N,1);

dT_non_lin(1) = heat_rate;

% forward differences in gradient
dT_non_lin(N1+N2+1:N-1) = ...
    (-lambda(N1+N2+1:N-1) ./ (rho(N1+N2+1:N-1) .* c_p(N1+N2+1:N-1).^2) .* dc_p(N1+N2+1:N-1) ...
     -lambda(N1+N2+1:N-1) ./ (rho(N1+N2+1:N-1).^2 .* c_p(N1+N2+1:N-1)) .* drho(N1+N2+1:N-1)) ...
    .* ((T(N1+N2+2:N) - T(N1+N2+1:N-1)).^2 ...
    ./ dx(N1+N2+1:N-1).^2);

dT_non_lin(N) = 0;


%% Linear part

% linear part vector
dT_lin = lambda(1:N) ./ dx.^2 ./ (c_p(1:N) .* rho(1:N)) .* (J_lin_sparse * T(1:N));

%% put linear and non-linear part together

dT = dT_non_lin + dT_lin;


return