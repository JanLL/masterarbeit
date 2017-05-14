N = 500;
L = 5.;
lambda = 0.96;

T_10 = simulate_1d('N', N, 'L', L, 'lambda', lambda, 'heat_rate', 10.);
T_5 = simulate_1d('N', N, 'L', L, 'lambda', lambda, 'heat_rate', 5.);
T_1 = simulate_1d('N', N, 'L', L, 'lambda', lambda, 'heat_rate', 1.);

plot(T_10(:,1), T_10(:,1)-T_10(:,end), 'DisplayName', 'beta=10'); hold on
plot(T_5(:,1), (T_5(:,1)-T_5(:,end))*1., 'DisplayName', 'beta=5'); hold on
plot(T_1(:,1), (T_1(:,1)-T_1(:,end))*1., 'DisplayName', 'beta=1'); hold on

xlabel('T_{oven}');
ylabel('T_{oven} - T_N');
legend(gca, 'show', 'Location', 'northwest')
