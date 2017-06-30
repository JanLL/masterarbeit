% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 10;
T_end = 200;

heat_rate = 10.; % K/min

lambda_test_setup = [23*1, 35.6000, 0.9600];

% measurement data
dsc = DSC204_readFile('ExpDat_16-407-3_mitKorr_10Kmin_H.csv');

% TODO: sinnvolles Intervall automatisch waehlen ...
index_T_29 = find(dsc.data(:,1) > 29, 1);
T_ref_meas = 30:0.05:157.9;


% evaluation points of simulation, every 0.05K, interpolation of
% measurements
U_dsc = interp1(dsc.data(index_T_29:end,1), dsc.data(index_T_29:end,3), ...
                T_ref_meas, 'pchip');
U_dsc = U_dsc * dsc.mass; % reverse normalization with mass [uV/mg] -> [uv]
U_dsc = transpose(cat(1, T_ref_meas, U_dsc));


% Solve optimization problem min_p ||U_dsc - dU||_2^2
knots = [-10,0, 30, 50, 60, 70, 80, 90, 100, 110, 115, 122, 127, 132, 135:155, 160, 170, 200];
%coeffs = 0.05 .* [1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 5, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
coeffs = 0.05 * ones(1, length(knots)-4);

k_sap_fit = [0.    45.5   0.]; % values from saphire-fit
k_data_table = [3.67763861e-02   6.00028439e+01  -4.47793211e+01]; % values from data table fit

p_optim_start = cat(2, knots, coeffs, k_sap_fit);

c_p_sample = {'B-', [length(knots), length(coeffs)]};
common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup,};
p_sim = get_param_sim(common_args{:}, 'c_p_sample', c_p_sample);
p_sim = update_c_p(p_sim, p_optim_start);

% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_estimable(end-2:end) = false; % fix mapping dT -> dU
p_optim_estimable(1:length(knots)) = false; % fix knot positions

p_optim_fixed = p_optim_start(~p_optim_estimable);

figure(1); % dU plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();

c_p_meas = calc_cp();

compute_residuum_expl = @(p_optim) ...
    compute_residuum(p_optim, p_optim_estimable, p_optim_fixed, p_sim, ...
                     U_dsc, c_p_meas, ax1, ax2);

% compute_residuum_expl(p_optim_start(p_optim_estimable)); % test initial value
% return


% knot_bounds = [-inf, 30, 70, 100, 120, 125, 130, 135, 140, 145, 150, inf];
% lb = cat(2, knot_bounds(1:end-1), ones(1, length(coeffs)+1)*-inf);
% ub = cat(2, knot_bounds(2:end), ones(1, length(coeffs)+1)*inf);

%lb = zeros(size(coeffs));
%ub = ones(size(coeffs))*100.;

lb = [];
ub = [];

opt_options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', 'OutputFcn', @disp_aux);
p_optim = lsqnonlin(compute_residuum_expl, p_optim_start(p_optim_estimable), lb, ub, opt_options);


p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

% Plot measured and optimized c_p
% figure();
% T_domain = linspace(dsc.data(1,1),dsc.data(end,1),200);
% % plot(T_domain, p_sim(1).eval_c_p(T_domain, p_sim(1).get_param_c_p(p_optim_all)), ...
% %      'DisplayName', 'Optimization'); hold on
% 
% % 
% % plot(dsc.data(:,1), c_p_meas, 'DisplayName', 'Measurement');
% 
% legend('show', 'location', 'northwest');
% xlabel('T_{ref}');
% ylabel('c_p');


% % print('/home/argo/masterarbeit/simulationen-data/c_p_optimized_lambda-wikiX8','-dpng');







