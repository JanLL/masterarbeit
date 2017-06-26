tic;

T0 = 30.;
heat_rate = 10. / 60;
L = 25.;

c_p = 0.41;
rho = 8.9;
lambda = 23.;

a = lambda / (c_p * rho);

x = L;
n = 10;


T_ref_max = 200';
t0 = 1/heat_rate*(T_ref_max - T0);
      
F = @(t) analytical_sol(x,t,n,T0, heat_rate, a) - T_ref_max;

fsolve_options = optimoptions('fsolve','Display','none');
t0_max = fsolve(F, t0, fsolve_options);

t = 0:1:t0_max;
T_ref_ana = analytical_sol(x,t,n,T0, heat_rate, a);


T_oven = T0 + heat_rate*t';
dT = T_oven - T_ref_ana;

toc

tic;
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

                
T_ref_num = simulate_1d(p_sim(1).eval_c_p, p_sim(1).eval_dc_p, p_sim(2));
toc



plot(T_ref_ana, dT, '--'); hold on
plot(T_ref_num(:,end), T_ref_num(:,1) - T_ref_num(:,end), '--');




