N = 500;
L = 0.005;
lambda = 0.96;


T_5 = ode_test(N, L, lambda, 5.);
T_1 = ode_test(N, L, lambda, 1.);

plot(T_5(:,1), T_5(:,1)-T_5(:,2)); hold on
plot(T_1(:,1), 5. .* (T_1(:,1)-T_1(:,2))); hold on