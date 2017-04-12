function f1 = dae_system1d(t, y)

startup;
t = 1.;
y = linspace(1,10,N)';
ae = zeros(N-1, 1);
y = cat(1, y, ae);

% y_in(1:N) =^ T(0), ..., T(N-1)
% y_in(N+1:end) =^ zusaetzl. DAE Gleichungen

% N variables for T, N-1 variables for aux. dae variables
f = zeros(N+N-1,1);

%% Non-linear part
f(1) = heat_rate;
f(2:N-1) = lambda / rho * -1./c_p_formula(y(2:N-1)).^2 .* ...
           dc_p_formula(y(2:N-1)) .* (y(3:N) - y(2:N-1)).^2 / dx^2;
f(N) = 0;


%% Linear part Jacobi matrix
F_lin = zeros(N,N);

% main diagonal 
F_lin = F_lin + diag( -2 ./ c_p_formula(y(1:N)));
F_lin(1,1) = 0.;
F_lin(N,N) = F_lin(N,N) / 2.;

% upper first diagonal
F_lin = F_lin + diag( 1 ./ c_p_formula(y(2:N)),1);
F_lin(1,2) = 0;

% lower first diagonal
F_lin = F_lin + diag( 1 ./ c_p_formula(y(2:N)),-1);


%% put linear and non-linear part together

f(1:N) = f(1:N) + lambda / (rho * dx^2) * F_lin*y(1:N);




return