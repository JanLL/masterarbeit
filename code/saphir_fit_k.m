% Measurement data

dsc_sap = DSC204_readFile('Sap-Kurve_10Kmin_H_Segment_7.csv');

index_T_start = find(dsc_sap.data(:,1) > 29, 1, 'first');
index_T_end = find(dsc_sap.data(:,1) < 157, 1, 'last');

U_dsc = [dsc_sap.data(index_T_start:index_T_end,1), dsc_sap.data(index_T_start:index_T_end,3)];
U_dsc(:,2) = U_dsc(:,2) * dsc_sap.mass; % reverse normalization with mass [uV/mg] -> [uv]


% Simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 10;
T_end = 200;

heat_rate = 10.; % K/min

% fit for different heat_rates assuming signal U_dsc grows linear with
% factor of heat_rate
U_dsc(:,2) = U_dsc(:,2) * heat_rate/10.;


lambda_test_setup = [23*1, 35.6000, 41.9]; % Wiki value: 41.9 
% lambda of saphire questionable

% coefficients of polynome to compute c_p of saphir in J/(g*K)
c_p_sap_coeffs = [  0.34407  ...
                   -0.37623  ...
                   -0.47824  ...
                    0.54579  ...
                    0.15393  ...
                   -0.10023  ...
                   -0.23778  ...
                    0.26410  ...
                   -0.21704  ...
                    0.23260  ...
                    1.12705 ];
      
                
c_p_sap = @(T, p) polyval(p, (T - 376.85) / 550);
dc_p_sap = @(T, p) polyval(cat(2, [0], p(1:end-1)), (T - 376.85) / 550);

c_p_sample = {c_p_sap, dc_p_sap, length(c_p_sap_coeffs)};
common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup ...
               'c_p_sample', c_p_sample};
p_sim = get_param_sim(common_args{:});
p_sim = update_c_p(p_sim, c_p_sap_coeffs);

k = [0., 1., 0.];
p_optim_start = cat(2, c_p_sap_coeffs, k);


% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_estimable(1:length(c_p_sap_coeffs)) = false;
p_optim_estimable(length(c_p_sap_coeffs) + 1) = false;
p_optim_estimable(length(c_p_sap_coeffs) + 3) = false;
p_optim_fixed = p_optim_start(~p_optim_estimable);

figure(1); % dU plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();

% pseudo measurement values of c_p to plot, equal to the ones used in
% simulation because we just optimize k here.
T_meas = (30:0.1:160)';
c_p_meas = [T_meas, p_sim(1).eval_c_p(T_meas)];

compute_residuum_expl = @(p_optim) ...
    compute_residuum(p_optim, p_optim_estimable, p_optim_fixed, p_sim, ...
                     U_dsc, c_p_meas, ax1, ax2);

% compute_residuum_expl(p_optim_start(p_optim_estimable)); % test initial value
% return


lb = [];
ub = [];

opt_options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', 'OptimalityTolerance', 1e-6);
[p_optim,~,~,~,optim_output] = lsqnonlin(compute_residuum_expl, p_optim_start(p_optim_estimable), lb, ub, opt_options);


% Comparison simulation fit and fit of data table from type_E_sensor.pdf
% For the latter see spielwiese.ipynb

data_table_fit_coeffs = [  5.09505329e-02   5.85548916e+01  -2.86468463e-01];

figure(3)
gca();
%cla;

T_domain_k = 0:0.1:30;

p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim;
p_optim_all(~p_optim_estimable) = p_optim_fixed;



plot(T_domain_k, polyval(p_sim.get_param_k(p_optim_all), T_domain_k), 'DisplayName', 'Optimization'); hold on
plot(T_domain_k, polyval(data_table_fit_coeffs, T_domain_k), 'DisplayName', 'Data table'); hold on
legend('show', 'location', 'northoutside')
xlabel('dT [K]')
ylabel('dU [uV]')
title('dU(dT) = k_0 + k_1 * dT + k_2 * dT^2')


p_optim



