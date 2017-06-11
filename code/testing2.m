L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 10;
T_end = 200;

lambda_test_setup = [23*1 35.6000 0.9600];


common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'lambda_test_setup', lambda_test_setup};
sim_options = get_param_sim(common_args{:}, 'heat_rate', 10.);


c_p_params_delta = [144.0009 - 15., ...
                    4.1036 * 5., ...
                    0.0039 + 0.1, ...
                    1.4217 * 0., ...
                    0.0078, ...
                    1.5325];

eval_c_p_expl = @(T) sim_options(1).eval_c_p(T, c_p_params_delta);
eval_dc_p_expl = @(T) sim_options(1).eval_dc_p(T, c_p_params_delta);


T_pcm_10 = simulate_1d(eval_c_p_expl, eval_dc_p_expl, sim_options(1));
T_ref_10 = simulate_1d(eval_c_p_expl, eval_dc_p_expl, sim_options(2));

N1 = sim_options(1).N1;

fig1 = figure(1);
dT10 = T_ref_10(:,N1) - T_pcm_10(:,N1);
plot(T_ref_10(:,N1), dT10, 'DisplayName', 'beta=10'); hold on