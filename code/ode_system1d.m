function dT = ode_system1d(t, T, N1, N2, N3, dx, heat_rate, ...
                           eval_c_p, eval_dc_p)
% [dT] = ode_system1d(t, T, N, dx, heat_rate, lambda, J_lin_sparse)
% 
% Computes the right hand side of the 1D differential heat equation for
% density and specific heat capacity is temperature dependent and thermal
% conductivity is constant.
%
% INPUT:     t --> time
%            T --> temperature in degree Celsius
%           N1 --> number of spatial discretization lattice points (Constantan)
%           N2 --> number of spatial discretization lattice points (Crucible)
%           N3 --> number of spatial discretization lattice points (PCM)
%           dx --> length [mm] of one spatial lattice point.
%    heat_rate --> rate [K/s] the temperature of the oven is increasing.
%     eval_c_p --> fhandle to evaluate specific heat capacity.
%    eval_dc_p --> fhandle to evaluate derivative of specific heat capacity
%                  w.r.t. temperature.
%
% OUTPUT:   dT --> right hand side of the 1D differential heat equation
%                  \nabla \left[ \frac{\lambda}{\rho c_p} \nabla T \right]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de

% initial definitions
N = N1+N2+N3;

% Check if J_lin_sparse was built, if not -> build it
persistent J_lin_sparse J_setup;

if isempty(J_lin_sparse) || isempty(J_setup) || any(J_setup ~= [N1, N2, N3])
    J_lin_sparse = build_linear_matrix(N);
    J_setup = [N1, N2, N3];
end


% pre-compute c_p and rho + derivatives because we need these later
% multiple times.
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
rho(N1+N2+1:end) = rho_formula(T(N1+N2+1:end)) .* 0.1;

drho = zeros(N, 1);
drho(N1+N2+1:end) = drho_formula(T(N1+N2+1:end)) .* 0.1;


lambda = ones(N, 1) * 0.96; % [mJ/mg*K], PCM, src: Robert: PCM_lambda.m
lambda(1:N1) = 23. * 8.; % [mW/(mm*K)], Constantan, src: Wikipedia
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