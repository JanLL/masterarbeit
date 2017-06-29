% values from properties_of_saphire.pdf

lambda_meas = [46.06, 25.12, 12.56]; % [mW/(mm K)]
temp_meas = [273, 373, 673]; % K

temp_meas = temp_meas - 273;  % K -> degC

calc_residuum = @(p) lambda_meas - (p(1) + p(2)./(temp_meas - p(3)));


p0 = [1,1,1];
opt_options = optimoptions('lsqnonlin', 'Display', 'iter-detailed');
lb = [];
ub = [];

p_optim = lsqnonlin(calc_residuum, p0, lb, ub, opt_options);


lambda_sap = @(T) p_optim(1) + p_optim(2)./(T - p_optim(3));
dlambda_sap = @(T) -p_optim(2)./(T-p_optim(3)).^2;

plot(temp_meas, lambda_meas, '.'); hold on

T_domain = -10:0.1:420;
plot(T_domain, lambda_sap(T_domain)); hold on
plot(T_domain, dlambda_sap(T_domain))


