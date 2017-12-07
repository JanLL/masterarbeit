function [residuum] = optimize_L1_residuum(optim_params, heat_rate, a_Const, T_0, dsc_measurements)
%

heat_rate_s = heat_rate / 60.;

T_ref_dsc = dsc_measurements.data(500:end,1);

nmp = length(T_ref_dsc);

L1 = optim_params(1);

n = 100;


F = @(t) analytical_sol(L1,t,n,T_0, heat_rate_s, a_Const);

dt = 1; % [s]

%residuum = F(t_0 + (0:1:(nmp-1))*dt) - T_ref_dsc;

meas_times = zeros(nmp,1);
for i=1:length(T_ref_dsc(:,1))

    F = @(t) analytical_sol(L1,t,n,T_0, heat_rate_s, a_Const) - T_ref_dsc(i);

    t_guess = (T_ref_dsc(i) - T_0)/heat_rate_s;
    %t_guess = 1.;
    
    fsolve_options = optimoptions('fsolve','Display','none');
    meas_times(i,1) = fsolve(F, t_guess, fsolve_options);
end

residuum = diff(meas_times) - 1.;


end

