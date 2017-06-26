T0 = 30.;
beta = 10. / 60;
L = 25.;

c_p = 0.41;
rho = 8.9;
lambda = 23.;

a = lambda / (c_p * rho);

x = 0:0.05:25;
t_end = 1/beta*(200 - T0);
t = 0:1:t_end;


% Analytical Solution
sol20 = analytical_sol(x,t,20,T0, beta, a);



% Numerical Solution
% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 200;

T_0 = 30;
T_end = 200;

heat_rate = 10.; % K/min

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate};
p_sim = get_param_sim(common_args{:});


c_p_params = [144.0009 - 15., ...
                    4.1036 * 5., ...
                    0.0039 + 0.1, ...
                    1.4217 * 0., ...
                    0.0078, ...
                    1.5325];

p_sim = update_c_p(p_sim, c_p_params);

                
T_ref = simulate_1d(p_sim(1).eval_c_p, p_sim(1).eval_dc_p, p_sim(2));


% Plot analytical and numerical result
plot(sol20(:,end), sol20(:,1) - sol20(:,end), '--'); hold on
plot(T_ref(:,end), T_ref(:,1) - T_ref(:,end), '--');
xlabel('T_{ref} [degC]');
ylabel('\Delta T');





