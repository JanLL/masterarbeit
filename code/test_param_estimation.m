% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 10;
T_end = 200;

heat_rate = 10.; % K/min

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate};
p_sim = get_param_sim(common_args{:});


% measurement data
dsc = DSC204_readFile('ExpDat_16-407-3_mitKorr_10Kmin_H.csv');

index_T_29 = find(dsc.data(:,1) > 29, 1);
T_ref_meas = 30:0.05:157.9;

% evaluation points of simulation, every 0.05K, interpolation of
U_dsc = interp1(dsc.data(index_T_29:end,1), dsc.data(index_T_29:end,3), ...
                T_ref_meas, 'pchip');

% Solve optimization problem min_p ||U_dsc - dU||_2^2

p_optim_start = [144.0009 - 15., ...
                 4.1036 * 5., ...
                 0.0039 + 0.005, ...
                 1.4217 * 1., ...
                 0.0078 + 0.01, ...
                 1.5325 + 4., ...
                 0.3];

compute_residuum_expl = @(p_optim) compute_residuum(p_optim, p_sim, U_dsc, T_ref_meas);

%compute_residuum_expl(p_optim_start);

lsqnonlin_options = {'Display', 'iter-detailed'};
p_optim = lsqnonlin(compute_residuum_expl, p_optim_start);


% Plot measured and optimized c_p
figure(2);
c_p_meas = calc_cp();
plot(dsc.data(:,1), c_p_meas, 'DisplayName', 'Measurement'); hold on

T_domain = linspace(dsc.data(1,1),dsc.data(end,1),200);
plot(T_domain, c_p_formula(T_domain, get_param_c_p(p_optim)), ...
     'DisplayName', 'Optimization');
legend('show', 'location', 'northwest');









