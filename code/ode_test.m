startup;


T0 = 30 .* ones(N,1);

t0 = 0.;
tf = 200.;
t = linspace(t0, tf, (tf-t0)*10000)';

opts = odeset('reltol', 1.0d-6, 'abstol', 1.0d-16);
sol = ode15s(@ode_system1d, t, T0, opts);

T = deval(sol, t)';

plot(t, T(:,1), 'g--', 'DisplayName', 'T_0'); hold on
plot(t, T(:,2), 'r--', 'DisplayName', 'T_1'); hold on
plot(t, T(:,3), 'y--', 'DisplayName', 'T_2'); hold on
plot(t, T(:,end), 'b--', 'DisplayName', 'T_N'); hold on


% xlabel('Time t');
% ylabel('Temp T');
% legend(gca, 'show', 'Location', 'northwest')

