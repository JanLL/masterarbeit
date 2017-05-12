function dT = ode_system1d(t, T, N, dx, heat_rate, lambda, J_lin_sparse)
% [dT] = ode_system1d(t, T, N, dx, heat_rate, lambda, J_lin_sparse)
% 
% Computes the right hand side of the 1D differential heat equation for
% density and specific heat capacity is temperature dependent and thermal
% conductivity is constant.
%
% INPUT:     t --> time
%            T --> temperature in degree Celsius.
%            N --> number of spatial discretization lattice points.
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


c_p = c_p_formula(T);
dc_p = dc_p_formula(T);

rho = rho_formula(T);
drho = drho_formula(T);


%% Non-linear part
dT_non_lin = zeros(N,1);

dT_non_lin(1) = heat_rate;


% backward differences in gradient
%dT_non_lin(2:N-1) = lambda / rho * -1./c_p(2:N-1).^2 .* ...
%           dc_p(2:N-1) .* (y(2:N-1) - y(1:N-2)).^2 / dx^2;


% forward differences in gradient
dT_non_lin(2:N-1) = ...
    (-lambda ./ (rho(2:N-1) .* c_p(2:N-1).^2) .* dc_p(2:N-1) ...
     -lambda ./ (rho(2:N-1).^2 .* c_p(2:N-1)) .* drho(2:N-1)) ...
    .* ((T(3:N) - T(2:N-1)).^2 / dx^2);

       
% central differences in gradient
%dT_non_lin(2:N-1) = lambda / rho * -1./c_p_formula(y(2:N-1)).^2 .* ...
%           dc_p_formula(y(2:N-1)) .* (y(3:N) - y(1:N-2)).^2 / (4*dx^2);


dT_non_lin(N) = 0;


%% Linear part

% linear part vector
dT_lin = lambda / dx^2 ./ (c_p(1:N) .* rho(1:N)) .* (J_lin_sparse * T(1:N));

%% put linear and non-linear part together

dT = dT_non_lin + dT_lin;


return