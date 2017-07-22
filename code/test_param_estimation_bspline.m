% measurement data
dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
dsc = DSC204_readFile(dsc_filename);

c_p_meas = calc_cp(dsc_filename);

% TODO: sinnvolles Intervall automatisch waehlen ... wobei das hier fuer
% alle Messungen bisher ganz gut war
index_T_dsc = [find(dsc.data(:,1) > 29, 1, 'first'), ...
               find(dsc.data(:,1) < 157.9, 1, 'last')];
U_dsc = [dsc.data(index_T_dsc(1):index_T_dsc(2),1), dsc.data(index_T_dsc(1):index_T_dsc(2),3)];

revMassNorm = true;  % reverse normalization with mass [uV/mg] -> [uv]
if revMassNorm
    U_dsc(:,2) = U_dsc(:,2) * dsc.mass;
end


% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;

N1 = 100;
N3 = 50;

T_0 = 10;
T_end = 200;

%heat_rate = 10.; % K/min
heat_rate = dsc.Tinfo.Tstep; % same heat rate as in measurement

lambda_test_setup = [23*1, 35.6000, 0.9600];

optim_solverName = 'lsqnonlin';

% Solve optimization problem min_p ||U_dsc - dU||_2^2
%knots = [-10,0, 30, 50, 60, 70, 80, 90, 100, 110, 115, 122, 127, 132, ...
%         135:1:160, 165, 170, 200];
knots = [-100, -50, 0, 10, 110, 115, 120, 125, 130, 132, 135:1:150, 155, 160, 165, 170, 220];
     
%coeffs = (0.05 .* [1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 5, 15, 15, 15, 15, 15, ...
%                  15, 15, 15, 15, 15, 15, 15, 15, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]);
coeffs = (0.02 .* [1, 3, 3, 3, 3, 3, 5, 15, 15, 15, 15, 15, ...
                  15, 15, 15, 15, 15, 15, 15, 15, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]);

% sqrt (later we square) to avoid negative coeffs when using optimizer without bounds

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
common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N1', N1, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup,};
p_sim = get_param_sim(common_args{:}, 'c_p_sample', c_p_sample);
p_sim = update_c_p(p_sim, p_optim_start);

% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_estimable(end-2:end) = false; % fix mapping dT -> dU
p_optim_estimable(1:length(knots)) = false; % fix knot positions
%p_optim_estimable(15:40) = true; % set knots at pos 135:160 free

n_free_knots = length(15:40);

p_optim_fixed = p_optim_start(~p_optim_estimable);

figure(1); % dU plot
ax1 = gca();
figure(2); % c_p plot
ax2 = gca();

if strcmp(optim_solverName, 'lsqnonlin')
    residuum_scalar_output = false;
elseif strcmp(optim_solverName, 'fminsearch') || ...
       strcmp(optim_solverName, 'fmincon')
    residuum_scalar_output = true;
end

compute_residuum_expl = @(p_optim) ...
    compute_residuum(p_optim, p_optim_estimable, p_optim_fixed, p_sim, ...
                     U_dsc, c_p_meas, residuum_scalar_output, ax1, ax2);

% TEST INITIAL VALUE
% compute_residuum_expl(p_optim_start(p_optim_estimable));
% return


if strcmp(optim_solverName, 'lsqnonlin')
    
    % Jacobian pattern generation
    Jac_pat = zeros(size(U_dsc,1),length(coeffs));

    for i=1:(length(coeffs))

        index_0 = find(U_dsc(:,1) > knots(max(1,i-20)),1,'first');
        if isempty(index_0)
            index_0 = 1;
        end

        index_1 = find(U_dsc(:,1) > knots(min(length(knots),i+4)),1,'first');
        if isempty(index_1)
            Jac_pat(index_0:end,i:end) = 1;
        end

        Jac_pat(index_0:index_1,i) = 1;    

    end
    
    
    lb = zeros(size(coeffs));
    ub = ones(size(coeffs))*100.;
    
    optim_con = {lb, ub};

    opt_options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
                               'OutputFcn', @disp_aux, ...
                               'StepTolerance', 1e-5, ...
                               'JacobPattern', Jac_pat);
    [p_optim,~,~,~,optim_output,~,jacobian] = lsqnonlin(...
        compute_residuum_expl, p_optim_start(p_optim_estimable), optim_con{:}, opt_options);

elseif strcmp(optim_solverName, 'fminsearch')
    optim_con = {};
    
    opt_options = optimset('Display', 'iter-detailed');
    [p_optim,~,~,optim_output] = fminsearch(compute_residuum_expl, p_optim_start(p_optim_estimable), opt_options);

elseif strcmp(optim_solverName, 'fmincon')
    
    n_knots = length(knots);
    n_coeffs = length(coeffs);

    A_columns = ones(n_free_knots-1, 2);

    % coeffs
    A_columns(:,2) = -1;
    
    diagonals = [0,1];
    A_sparse = spdiags(A_columns, diagonals, n_free_knots-1, n_free_knots+n_coeffs);

    A = full(A_sparse);
    b = zeros(size(A,1),1);

    lb = [135 * ones(1,n_free_knots), zeros(1,n_coeffs)];
    ub = [160 * ones(1,n_free_knots), ones(1,n_coeffs)*100.];
    
    Aeq = [];
    beq = [];
    %lb = zeros(1,n_coeffs);
    %ub = ones(1,n_coeffs)*100.;
    nonlincon = @nonlcon_empty;
    
    optim_con = {A, b, Aeq, beq, lb, ub, nonlincon};

    opt_options = optimoptions('fmincon', 'Display', 'iter-detailed', 'OutputFcn', @disp_aux);
    [p_optim,~,~,optim_output] = fmincon(compute_residuum_expl, p_optim_start(p_optim_estimable), ...
                      optim_con{:}, opt_options);

else
    error('Choose optim_solverName from [''lsqnonlin'', ''fminsearch'', ''fmincon'']!');
end
    

% update all (free and fixed) optimization parameters with optimized values
p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

% Plot final results
compute_residuum_expl(p_optim_all(p_optim_estimable));

p_sim = update_c_p(p_sim, p_optim_all);

save_path = '/home/argo/masterarbeit/fits_data/';
fit_data = save_fit(save_path, dsc, index_T_dsc, revMassNorm, ...
    p_sim, optim_solverName, opt_options, p_optim_start, p_optim_estimable, ...
    optim_con, p_optim_all, optim_output);

save('/home/argo/masterarbeit/fits_data/jacobian.mat', 'jacobian')
