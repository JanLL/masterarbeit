% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 10;
T_end = 200;

heat_rate = 10.; % K/min

lambda_test_setup = [23*4 35.6000 0.9600];

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup,};
p_sim = get_param_sim(common_args{:});


% measurement data
dsc = DSC204_readFile('ExpDat_16-407-3_mitKorr_10Kmin_H.csv');

% TODO: sinnvolles Intervall automatisch waehlen ...
index_T_29 = find(dsc.data(:,1) > 29, 1);
T_ref_meas = 30:0.05:157.9;

% evaluation points of simulation, every 0.05K, interpolation of
% measurements
U_dsc = interp1(dsc.data(index_T_29:end,1), dsc.data(index_T_29:end,3), ...
                T_ref_meas, 'pchip');

% Solve optimization problem min_p ||U_dsc - dU||_2^2
knots = [-0, 50., 90., 110., 122, 127, 132, 137, 142, 147, 200];
coeffs = [10, 15, 20., 100., 1.5, 1.5, 1.5];
k = 0.3;
p_optim_start = cat(2, knots, coeffs, k);

c_p_sample = {'B-', [length(knots), length(coeffs)]};
common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup,};
p_sim = get_param_sim(common_args{:}, 'c_p_sample', c_p_sample);


% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_estimable(1:length(knots)) = false;
p_optim_fixed = p_optim_start(~p_optim_estimable);

compute_residuum_expl = @(p_optim) ...
    compute_residuum(p_optim, p_optim_estimable, p_optim_fixed, p_sim, U_dsc, T_ref_meas);


% compute_residuum_expl(p_optim_start(p_optim_estimable)); % test initial value
% return


knot_bounds = [-inf, 30, 70, 100, 120, 125, 130, 135, 140, 145, 150, inf];

%lb = cat(2, knot_bounds(1:end-1), ones(1, length(coeffs)+1)*-inf);
%ub = cat(2, knot_bounds(2:end), ones(1, length(coeffs)+1)*inf);
lb = [];
ub = [];


opt_options = optimoptions('lsqnonlin', 'Display', 'iter-detailed');

p_optim = lsqnonlin(compute_residuum_expl, p_optim_start(p_optim_estimable), lb, ub, opt_options);


p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim;
p_optim_all(~p_optim_estimable) = p_optim_fixed;


% Plot measured and optimized c_p
figure();
T_domain = linspace(dsc.data(1,1),dsc.data(end,1),200);
plot(T_domain, p_sim(1).eval_c_p(T_domain, p_sim(1).get_param_c_p(p_optim_all)), ...
     'DisplayName', 'Optimization'); hold on

[T_meas, c_p_meas] = calc_cp();
plot(dsc.data(:,1), c_p_meas, 'DisplayName', 'Measurement');

legend('show', 'location', 'northwest');
xlabel('T_{ref}');
ylabel('c_p');


% % print('/home/argo/masterarbeit/simulationen-data/c_p_optimized_lambda-wikiX8','-dpng');







