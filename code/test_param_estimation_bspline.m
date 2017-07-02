% measurement data
dsc = DSC204_readFile('ExpDat_16-407-3_mitKorr_10Kmin_H.csv');

c_p_meas = calc_cp(); % TODO: calc_cp allgemein fuer beliebige Messdaten als input

% TODO: sinnvolles Intervall automatisch waehlen ...
index_T_dsc = [find(dsc.data(:,1) > 29, 1, 'first'), ...
               find(dsc.data(:,1) < 157.9, 1, 'last')];
U_dsc = [dsc.data(index_T_dsc(1):index_T_dsc(2),1), dsc.data(index_T_dsc(1):index_T_dsc(2),3)];

revMassNorm = false;  % reverse normalization with mass [uV/mg] -> [uv]
if revMassNorm
    U_dsc(:,2) = U_dsc(:,2) * dsc.mass;
end


% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 10;
T_end = 200;

heat_rate = 10.; % K/min

lambda_test_setup = [23*1, 35.6000, 0.9600];

optim_solverName = 'lsqnonlin';
% optim_solverName = 'fminsearch';
% optim_solverName = 'fmincon';
% Solve optimization problem min_p ||U_dsc - dU||_2^2
knots = [-10,0, 30, 50, 60, 70, 80, 90, 100, 110, 115, 122, 127, 132, 135:160, 165, 170, 200];
coeffs = (0.05 .* [1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 5, 15, 15, 15, 15, 15, ...
                  15, 15, 15, 15, 15, 15, 15, 15, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]);
% sqrt (later we square) to avoid negative coeffs when using optimizer without bounds
              
%coeffs = 0.05 * ones(1, length(knots)-4);

% length(knots)
% length(coeffs)
% return

if revMassNorm
    k_sap_fit = [0.    45.5951   0.]; % values from saphire-fit
else
    k_sap_fit = [0.    0.5381   0.]; % values from saphire-fit
end
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


compute_residuum_expl = @(p_optim) ...
    compute_residuum(p_optim, p_optim_estimable, p_optim_fixed, p_sim, ...
                     U_dsc, c_p_meas, ax1, ax2);

% TEST INITIAL VALUE
% compute_residuum_expl(p_optim_start(p_optim_estimable));
% return


% knot_bounds = [-inf, 30, 70, 100, 120, 125, 130, 135, 140, 145, 150, inf];
% lb = cat(2, knot_bounds(1:end-1), ones(1, length(coeffs)+1)*-inf);
% ub = cat(2, knot_bounds(2:end), ones(1, length(coeffs)+1)*inf);

if strcmp(optim_solverName, 'lsqnonlin')
    lb = zeros(size(coeffs));
    ub = ones(size(coeffs))*100.;
    
    optim_con = {lb, ub};

    opt_options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', 'OutputFcn', @disp_aux);
    [p_optim,~,~,~,optim_output] = lsqnonlin(...
        compute_residuum_expl, p_optim_start(p_optim_estimable), optim_con{:}, opt_options);

elseif strcmp(optim_solverName, 'fminsearch')
    optim_con = {};
    
    opt_options = optimset('Display', 'iter-detailed');
    p_optim = fminsearch(compute_residuum_expl, p_optim_start(p_optim_estimable), opt_options);

elseif strcmp(optim_solverName, 'fmincon')
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = zeros(size(coeffs));
    ub = ones(size(coeffs))*100.;
    nonlincon = @nonlcon_empty;
    
    optim_con = {A, b, Aeq, beq, lb, ub, nonlincon};

    opt_options = optimoptions('fmincon', 'Display', 'iter-detailed', 'OutputFcn', @disp_aux);
    p_optim = fmincon(compute_residuum_expl, p_optim_start(p_optim_estimable), ...
                      optim_con{:}, opt_options);

else
    error('Choose optim_solverName from [lsqnonlin, fminsearch, fmincon]!');
end
    

% update all (free and fixed) optimization parameters with optimized values
p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

% Plot final results
compute_residuum_expl(p_optim_all(p_optim_estimable));

p_sim = update_c_p(p_sim, p_optim_all);

save_path = '/home/argo/masterarbeit/simulationen-data/test/';
save_fit(save_path, dsc, index_T_dsc, revMassNorm, ...
    p_sim, optim_solverName, opt_options, p_optim_start, p_optim_estimable, ...
    optim_con, p_optim_all, optim_output);


