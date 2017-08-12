% measurement data
dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
dsc = DSC204_readFile(dsc_filename);

m_pcm = dsc.mass;
c_p_meas = calc_cp(dsc_filename);

% TODO: sinnvolles Intervall automatisch waehlen ... wobei das hier fuer
% alle Messungen bisher ganz gut war
index_T_dsc = [find(dsc.data(:,1) > 29, 1, 'first'), ...
               find(dsc.data(:,1) < 157.9, 1, 'last')];
q_dsc = [dsc.data(index_T_dsc(1):index_T_dsc(2),1), ...
         dsc.data(index_T_dsc(1):index_T_dsc(2),3) ...
         ./ dsc.data(index_T_dsc(1):index_T_dsc(2),4)];

revMassNorm = true;  % reverse normalization with mass [uV/mg] -> [uv]
if revMassNorm
    q_dsc(:,2) = q_dsc(:,2) * m_pcm;
end


% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 10;
T_end = 200;

%heat_rate = 10.; % K/min
heat_rate = dsc.Tinfo.Tstep; % same heat rate as in measurement

lambda_test_setup = [23*1, 35.6000, 0.9600];

optim_solverName = 'lsqnonlin';
optim_type = 'heat_flux_1';
optim_type_int = str2double(optim_type(end));

% c_p parametrization of sample
nrb_order = 4; % nrb_order = 4 equates to C^2

cntrl_pts_x = [0, 30, 60, 90, 100:2:150, 160, 180, 200]; %24
cntrl_pts_y = [repmat(2.,1,length(cntrl_pts_x))];
cntrl_pts = [cntrl_pts_x; cntrl_pts_y];
num_cntrl_pts = size(cntrl_pts,2);


% equidistant knots
knots = [zeros(1,nrb_order), ...
         (1:num_cntrl_pts-nrb_order)/(num_cntrl_pts-nrb_order+1), ...
         ones(1,nrb_order)];

c_p_sample = {'NURBS', [num_cntrl_pts, length(knots)]};
common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup, ...
               'c_p_sample', c_p_sample};
     
p_sim = get_param_sim(common_args{:}, 'c_p_sample', c_p_sample);

p_optim_start = [cntrl_pts(1,:), cntrl_pts(2,:), knots];

p_sim = update_c_p(p_sim, p_optim_start);

% choose free(true)/fixed(false) parameters to optimize
p_optim_estimable = true(length(p_optim_start), 1);
p_optim_estimable(1:num_cntrl_pts) = false; % fix x-position of control pts
p_optim_estimable(2*num_cntrl_pts+1:2*num_cntrl_pts+length(knots)) = false; % fix knot vector
%p_optim_estimable(15:40) = true; % set knots at pos 135:160 free

p_optim_fixed = p_optim_start(~p_optim_estimable);

figure(1); % q_pcm_in plot
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
    compute_residuum_q(p_optim, optim_type_int, p_optim_estimable, p_optim_fixed, p_sim, ...
                       q_dsc, c_p_meas, m_pcm, residuum_scalar_output, ax1, ax2);


% TEST INITIAL VALUE
% compute_residuum_expl(p_optim_start(p_optim_estimable));
% return

% LOAD FIT RESULT
% fit_data = load('fit_data.mat');
% p_optim_all = fit_data.optimization.param_end;
% 
% p_sim = update_c_p(fit_data.simulation, p_optim_all);
% 
% % TODO: compute_residuum_expl neu aufbauen mit daten aus fit_data um
% % q_pcm_in und c_p neu zu berechnen und plotten...
% 
% T_pcm = simulate_1d(fit_data.simulation);
% 
% max(T_pcm(:,end-50) - T_pcm(:,end))
% 
% figure(3)
% image(T_pcm(:,end-50:end), 'CDataMapping', 'scaled')
% colorbar;
% colormap hsv;
% title('Temp. distribution in PCM');
% ylabel('time');
% xlabel('space [grid #]');
% grid;
% 
% return


if strcmp(optim_solverName, 'lsqnonlin')
    lb = zeros(1,num_cntrl_pts);
    ub = ones(1,num_cntrl_pts)*200.;
    optim_con = {lb, ub};

    opt_options = optimoptions('lsqnonlin', ...
                               'Display', 'iter-detailed', ...
                               'OutputFcn', @disp_aux);%, ...
                               %'SpecifyObjectiveGradient', true);
    [p_optim,~,~,~,optim_output,~,jac_output] = lsqnonlin(...
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
save_commong_args = {save_path, dsc, index_T_dsc, revMassNorm, p_sim, ...
    optim_type, optim_solverName, opt_options, p_optim_start, ...
    p_optim_estimable, optim_con, p_optim_all, optim_output};

if strcmp(optim_solverName, 'lsqnonlin')
    fit_data = save_fit(save_commong_args{:}, jac_output);
elseif strcmp(optim_solverName, 'fminsearch') || ...
       strcmp(optim_solverName, 'fmincon')
    fit_data = save_fit(save_commong_args{:});
end
