startup;


T0 = 90 .* ones(N,1);

t0 = 0.;
tf = 25000.;

t_domain = linspace(t0, tf, (tf-t0));

[t,T] = ode23t(@(t,y) ode_system1d(t,y), [t0 tf], T0);

plot(t, T(:,1), 'g--', 'DisplayName', 'T_0'); hold on
plot(t, T(:,2), 'r--', 'DisplayName', 'T_1'); hold on
plot(t, T(:,3), 'y--', 'DisplayName', 'T_2'); hold on

plot(t, T(:,end), 'b--', 'DisplayName', 'T_N'); hold on

xlabel('Time t');
ylabel('Temp T');
legend(gca, 'show', 'Location', 'northwest')