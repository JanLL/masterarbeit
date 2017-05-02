startup;


T0 = 80 .* ones(N,1);

t0 = 0.;
tf = 50.;
t = linspace(t0, tf, (tf-t0)*10000)';


%Jacobi matrix sparse variant
J_lin_columns = ones(N, 3);

% main diagonal
J_lin_columns(2:end-1, 1) = -2.;
J_lin_columns(1, 1) = 0;
J_lin_columns(N, 1) = -1.;

% upper first diagonal
J_lin_columns(1, 2) = 0.;  % one element less than in main diagonal

% lower first diagonal
J_lin_columns(1:2, 2) = 0.;  % one element less again and one zero entry

% build actual sparse matrix for linear part
diagonals = [0, 1, -1];
J_lin_sparse = spdiags(J_lin_columns, diagonals, N, N);


opts = odeset('reltol', 1.0d-6, 'abstol', 1.0d-16);

tic;
sol = ode15s(@(t, y) ode_system1d(t, y, J_lin_sparse), t, T0, opts);
toc

T = deval(sol, t)';

plot(t, T(:,1) - T(:,1), 'g--', 'DisplayName', 'T_0'); hold on
plot(t, T(:,2) - T(:,1), 'r--', 'DisplayName', 'T_1'); hold on
plot(t, T(:,3) - T(:,1), 'y--', 'DisplayName', 'T_2'); hold on
%plot(t, T(:,end), 'b--', 'DisplayName', 'T_N'); hold on


% xlabel('Time t');
% ylabel('Temp T');
% legend(gca, 'show', 'Location', 'northwest')

