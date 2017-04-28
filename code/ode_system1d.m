function dT = ode_system1d(t, T)

startup;


%t = 1.;
%y = linspace(1,10,N)';

c_p = c_p_formula(T);
dc_p = dc_p_formula(T);

% c_p = 5 * ones(size(c_p));
% dc_p = 0 * ones(size(c_p));

dT = zeros(N,1);

%% Non-linear part
dT_non_lin = zeros(N,1);

dT_non_lin(1) = heat_rate;


% backward differences in gradient
%dT_non_lin(2:N-1) = lambda / rho * -1./c_p(2:N-1).^2 .* ...
%           dc_p(2:N-1) .* (y(2:N-1) - y(1:N-2)).^2 / dx^2;

% forward differences in gradient
dT_non_lin(2:N-1) = lambda / rho * -1./c_p(2:N-1).^2 .* ...
           dc_p(2:N-1) .* (T(3:N) - T(2:N-1)).^2 / dx^2;

% central differences in gradient
%dT_non_lin(2:N-1) = lambda / rho * -1./c_p_formula(y(2:N-1)).^2 .* ...
%           dc_p_formula(y(2:N-1)) .* (y(3:N) - y(1:N-2)).^2 / (4*dx^2);


dT_non_lin(N) = 0;


%% Linear part
%Jacobi matrix

% todo: sparse:
% spdiag, sparse(N,N)
% siehe spdiags(B,d,A)

J_lin = zeros(N,N);

% main diagonal 
J_lin = J_lin + diag(-2. * ones(N,1));
J_lin(1,1) = 0.;
J_lin(N,N) = -1.;

% upper first diagonal
J_lin = J_lin + diag(1. * ones(N-1,1),1);
J_lin(1,2) = 0;

% lower first diagonal
J_lin = J_lin + diag(1. * ones(N-1,1),-1);

% inverse c_p vector
inv_c_p = 1 ./ c_p_formula(T(1:N));

% linear part vector
dT_lin = lambda / (rho * dx^2) * inv_c_p .* (J_lin * T(1:N));

%% put linear and non-linear part together

dT = dT_non_lin + dT_lin;


return