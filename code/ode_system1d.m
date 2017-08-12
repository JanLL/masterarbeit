function dT = ode_system1d(t, T, p_sim, dx)
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
% c_p_test_setup --> 2x1 array: specific heat capacity [mJ/(mg*K] 
%                               of [Constantan, crucible]
% rho_test_setup --> 2x1 array: density [mg/mm^3] of [Constantan, crucible]
%         lambda --> 3x1 array: heat conductivity [mW/(mm*K)] 
%                               of [Constantan, crucible, PCM]
%     eval_c_p --> fhandle to evaluate specific heat capacity.
%    eval_dc_p --> fhandle to evaluate derivative of specific heat capacity
%                  w.r.t. temperature.
%
% OUTPUT:   dT --> right hand side of the 1D differential heat equation
%                  \nabla \left[ \frac{\lambda}{\rho c_p} \nabla T \right]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de

% initial definitions
N1 = p_sim.N1;
N2 = p_sim.N2;
N3 = p_sim.N3;
N = N1+N2+N3;

L1 = p_sim.L1;
L2 = p_sim.L2;
L3 = p_sim.L3;

heat_rate = p_sim.heat_rate / 60; % [K/min] -> [K/min]

% Check if J_lin_sparse was built, if not -> build it
persistent J_lin_sparse J_setup;

if isempty(J_lin_sparse) || isempty(J_setup) || ...
   any(J_setup ~= [N1, N2, N3, L1, L2, L3])
    J_lin_sparse = build_linear_matrix(N);
    
    alpha = L3/L1 * N1/N3;

    % Constantan - PCM transition with inhomogenous grid
    J_lin_sparse(N1,N1-1) = 2/(1+alpha);
    J_lin_sparse(N1,N1) = -2/alpha;
    J_lin_sparse(N1,N1+1) = 2/(alpha*(alpha+1));
    
    J_setup = [N1, N2, N3, L1, L2, L3];
end

% pre-compute c_p and rho + derivatives because we need these later
% multiple times.
c_p = ones(N, 1);
c_p(1:N1) = p_sim.c_p_test_setup(1); 
c_p(N1+1:N1+N2) = p_sim.c_p_test_setup(2); 
c_p(N1+N2+1:end) = p_sim.eval_c_p(T(N1+N2+1:end));

dc_p = zeros(N, 1);
dc_p(N1+N2+1:end) = p_sim.eval_dc_p(T(N1+N2+1:end));


rho = ones(N, 1);
rho(1:N1) = p_sim.rho_test_setup(1); 
rho(N1+1:N1+N2) = p_sim.rho_test_setup(2); 
rho(N1+N2+1:end) = rho_formula(T(N1+N2+1:end));% .* 0.1;

% factor 0.1 (a simple guess) to compensate different 
% masses(-> cross sections) of Constantan/PCM.

drho = zeros(N, 1);
drho(N1+N2+1:end) = drho_formula(T(N1+N2+1:end));% .* 0.1;

lambda = ones(N, 1);
lambda(1:N1) = p_sim.lambda_test_setup(1);
lambda(N1+1:N1+N2) = p_sim.lambda_test_setup(2);
lambda(N1+N2+1:end) = p_sim.lambda_test_setup(3);

%% Non-linear part
dT_non_lin = zeros(N,1);

dT_non_lin(1) = heat_rate;

% forward differences in gradient
dT_non_lin(N1+N2+1:N-1) = ...
    (-lambda(N1+N2+1:N-1) ./ (rho(N1+N2+1:N-1) .* c_p(N1+N2+1:N-1).^2) .* dc_p(N1+N2+1:N-1) ...
     -lambda(N1+N2+1:N-1) ./ (rho(N1+N2+1:N-1).^2 .* c_p(N1+N2+1:N-1)) .* drho(N1+N2+1:N-1)) ...
    .* ((T(N1+N2+2:N) - T(N1+N2+1:N-1)).^2 ...
    ./ dx(N1+N2+2:N).^2);

dT_non_lin(N) = 0;


%% Linear part

% linear part vector
dT_lin = lambda(1:N) ./ dx.^2 ./ (c_p(1:N) .* rho(1:N)) .* (J_lin_sparse * T(1:N));

%% put linear and non-linear part together

dT = dT_non_lin + dT_lin;


return