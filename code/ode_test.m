function [T] = ode_test(N, L, lambda, heat_rate)

if ~exist('N','var'), N=500; end
if ~exist('L','var'), L=0.005; end %[m]
if ~exist('lambda','var'), lambda=0.96; end % [W/(m * K)]
if ~exist('heat_rate','var'), heat_rate=5.; end % [K/min]

dx = L/N;

% integration initial values
T0 = 80 .* ones(N,1);

t0 = 0.;
tf = 150. / heat_rate;
t = linspace(t0, tf, (tf-t0)*100)';

% pre-compute constant sparse matrix from linear part
J_lin_sparse = build_linear_matrix(N);


opts = odeset('reltol', 1.0d-8, 'abstol', 1.0d-16);

tic;
sol = ode15s(@(t, y) ode_system1d(t, y, N, dx, heat_rate, lambda, J_lin_sparse), t, T0, opts);
toc


T = deval(sol, t)';


% figure(1);
% plot(t, T(:,1) - T(:,1), 'g', 'DisplayName', 'T_0'); hold on
% plot(t, T(:,2) - T(:,1), 'r', 'DisplayName', 'T_1');
% plot(t, T(:,3) - T(:,1), 'y', 'DisplayName', 'T_2');
% plot(t, T(:,end) - T(:,1), 'b', 'DisplayName', 'T_N');

%figure(2);
%plot(T(:,1), T(:,1)-T(:,2), 'DisplayName', 'T_0 - T_1'); hold on

% xlabel('Time t');
% ylabel('Temp T');
% legend(gca, 'show', 'Location', 'northwest')


end

