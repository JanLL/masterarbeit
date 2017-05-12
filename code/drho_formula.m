function [drho] = drho_formula(T)
% [rho] = drho_formula(T)
% 
% Calculates density for specified temperature.
%
% INPUT:    T --> temperature in [degree Celsius]
%
% OUTPUT: drho --> d(density)/dT in [mg/mm^3]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de



% formula and parameter from PCM_rho.m from Robert
p0 = 167.2182;
p1 = 129.9861;
p2 = 760.8218;
p3 = 0.078916;


drho = -p0*p3*(exp(p3*(T-p1))) ./ (1+exp(p3*(T-p1))).^2;
drho = drho * 1e-3;

return