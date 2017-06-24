% normal forward integration


L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 30;
T_end = 200;

heat_rate = 10.; % K/min

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate};
p_sim = get_param_sim(common_args{:});


c_p_params = [144.0009, ...
                    4.1036, ...
                    0.0039, ...
                    1.4217, ...
                    0.0078, ...
                    1.5325];

p_sim = update_c_p(p_sim, c_p_params);


T_pcm = simulate_1d(p_sim(1).eval_c_p, p_sim(1).eval_dc_p, p_sim(1));
T_ref = simulate_1d(p_sim(1).eval_c_p, p_sim(1).eval_dc_p, p_sim(2));


dT = T_ref(:,p_sim(1).N1) - T_pcm(:,p_sim(1).N1);

plot(T_ref(:,p_sim(1).N1), dT)

