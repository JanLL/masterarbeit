L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 10;
T_end = 200;

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end};
sim_options_10 = get_param_sim(common_args{:}, 'heat_rate', 10.);
sim_options_5 = get_param_sim(common_args{:}, 'heat_rate', 5.);
sim_options_1 = get_param_sim(common_args{:}, 'heat_rate', 1.);

% set type and parameters of c_p
c_p_params_delta = [144.0009 - 15., ...
                    4.1036 * 5., ...
                    0.0039 + 0.1, ...
                    1.4217 * 0., ...
                    0.0078, ...
                    1.5325];
c_p_params_meas_10K = [144.0009, ...
                        4.1036, ...
                        0.0039, ...
                        1.4217, ...
                        0.0078, ...
                        1.5325];
c_p_params_test = [144.0009 - 15., ...
                   4.1036 * 5., ...
                   0.0039 + 0.005, ...
                   1.4217 * 1., ...
                   0.0078 + 0.01, ...
                   1.5325 + 4];
                    

% compute derivative (symbolically) of c_p w.r.t. temperature T
syms T;
p = sym('p', [6 1]);
dc_p = matlabFunction(diff(c_p_formula(T, p), T), 'Vars', [T;p]);

% build new function handles of c_p and derivative with explicit parameter
% values
c_p_params = c_p_params_test;
eval_c_p = @(T)c_p_formula(T, c_p_params);
c_p_params_cell = num2cell(c_p_params);
eval_dc_p = @(T) dc_p(T,c_p_params_cell{1:end});

% TODO: same for rho ...

% PCM side 
T_pcm_10 = simulate_1d(eval_c_p, eval_dc_p, sim_options_10(1));
% T_pcm_5 = simulate_1d(eval_c_p, eval_dc_p, sim_options_5(1));
% T_pcm_1 = simulate_1d(eval_c_p, eval_dc_p, sim_options_1(1));


% Reference side
T_ref_10 = simulate_1d(eval_c_p, eval_dc_p, sim_options_10(2));
% T_ref_5 = simulate_1d(eval_c_p, eval_dc_p, sim_options_5(2));
% T_ref_1 = simulate_1d(eval_c_p, eval_dc_p, sim_options_1(2));


% take Delta T between reference and link constantan-pcm and plot against T_ref
N1 = sim_options_10.N1;

fig1 = figure(1);
dT10 = T_ref_10(:,N1) - T_pcm_10(:,N1);
plot(T_ref_10(:,N1), dT10, 'DisplayName', 'beta=10'); hold on

plot(T_ref_meas, 3.*U_dsc, 'DisplayName', 'Measurements')

legend('show', 'location', 'northwest');

% dT5 = T_ref_5(:,N1) - T_pcm_5(:,N1);
% plot(T_ref_5(:,N1), 2.*dT5, 'DisplayName', 'beta=5'); hold on
% 
% dT1 = T_ref_1(:,N1) - T_pcm_1(:,N1);
% plot(T_ref_1(:,N1), 10.*dT1, 'DisplayName', 'beta=1'); hold on



           