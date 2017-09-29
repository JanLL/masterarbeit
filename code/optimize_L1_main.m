

lambda_Const = 23.; 
rho_Const = 8.9;    
c_p_Const = 0.41;   

a_Const = lambda_Const / (rho_Const * c_p_Const);

T_0 = -50.; % [degC]



% Measurements
dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
dsc_measurement = DSC204_readFile(dsc_filename);

heat_rate = dsc_measurement.Tinfo.Tstep; % [K/min]


t_0 = 500;
T_0 = 10.;



optimize_L1_residuum_expl = @(optim_params) optimize_L1_residuum(optim_params, heat_rate, a_Const, T_0, dsc_measurement);


% start values
L1 = 10.;

p_optim_start = [L1];

% Start initial value
% sum(optimize_L1_residuum_expl(p_optim_start).^2)
% return



opt_options = optimoptions('lsqnonlin', ...
                           'Display', 'iter-detailed', ...
                           'MaxIter', 1000, ...
                           'StepTolerance', 1e-6, ...
                           'FunctionTolerance', 1e-10);

lb = [];
ub = [];
opt_constraints = {lb, ub};


[p_optim_end,~,~,~,optim_output] = lsqnonlin(optimize_L1_residuum_expl, p_optim_start, opt_constraints{:}, opt_options);
%optimize_L1_residuum_expl(p_optim_end)
p_optim_end

return


% Test measurement time points
L1 = p_optim_end(1);

heat_rate_s = heat_rate / 60.;

T_ref_dsc = dsc_measurement.data(500:end,1);
nmp = length(T_ref_dsc);


% Analytical solution for T_ref: Solve non-linear system to get time t 
% where T_ref (temp. at crucible) reaches T_ref_meas (time where heat flux
% measurement was done)
n = 100;
meas_times = zeros(nmp, 1);


for i=1:length(T_ref_dsc(:,1))

    F = @(t) analytical_sol(L1,t,n,T_0, heat_rate_s, a_Const) - T_ref_dsc(i);

    t_guess = (T_ref_dsc(i) - T_0)/heat_rate_s;
    %t_guess = 1.;
    
    fsolve_options = optimoptions('fsolve','Display','none');
    meas_times(i,1) = fsolve(F, t_guess, fsolve_options);
end



plot(diff(meas_times(1:end)))


meas_times(1:20)







