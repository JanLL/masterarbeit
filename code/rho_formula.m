function [rho] = rho_formula(T)

% formula from PCM_rho.m from Robert

% rho in [kg/m^3] ??
% T in [°C]

p0 = 167.2182;
p1 = 129.9861;
p2 = 760.8218;
p3 = 0.078916;


rho = p0./(1.0+exp(p3.*(T-p1))) + p2;

return