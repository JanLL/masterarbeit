N1 = 50;
L1 = 50.;
N2 = 450;
L2 = 5.;
lambda = 0.96;

T_end = 600;

T_10 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                   'lambda', lambda, 'heat_rate', 10., 'T_end', T_end);
T_5 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                  'lambda', lambda, 'heat_rate', 5., 'T_end', T_end);
T_1 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                  'lambda', lambda, 'heat_rate', 1., 'T_end', T_end);

plot(T_10(:,1), T_10(:,N1-1)-T_10(:,N1), 'DisplayName', 'beta=10'); hold on
plot(T_5(:,1), (T_5(:,N1-1)-T_5(:,N1))*1., 'DisplayName', 'beta=5'); hold on
plot(T_1(:,1), (T_1(:,N1-1)-T_1(:,N1))*1., 'DisplayName', 'beta=1'); hold on

xlabel('T_{oven}');
ylabel('T_{oven} - T_N');
legend(gca, 'show', 'Location', 'northwest')
