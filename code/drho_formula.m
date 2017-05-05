function [drho] = drho_formula(T)

% formula from PCM_rho.m from Robert

% rho in [kg/m^3] ??
% T in [°C]

p0 = 167.2182;
p1 = 129.9861;
p2 = 760.8218;
p3 = 0.078916;


drho = -p0*p3*(exp(p3*(T-p1))) ./ (1+exp(p3*(T-p1))).^2;