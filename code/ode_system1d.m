function dT = ode_system1d(t, T, N, dx, heat_rate, lambda, J_lin_sparse)


%t = 1.;
%y = linspace(1,10,N)';

c_p = c_p_formula(T);
dc_p = dc_p_formula(T);

rho = rho_formula(T);
drho = drho_formula(T);

% test case where rho is constant
%rho = ones(size(rho)) .* 800;
%drho = drho .* 0;


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
% inverse c_p vector
%inv_c_p = 1 ./ (c_p(1:N) .* rho(1:N));

% linear part vector
dT_lin = lambda / dx^2 ./ (c_p(1:N) .* rho(1:N)) .* (J_lin_sparse * T(1:N));

%% put linear and non-linear part together

dT = dT_non_lin + dT_lin;


return