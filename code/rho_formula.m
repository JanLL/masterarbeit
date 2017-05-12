function [rho] = rho_formula(T)
% [rho] = rho_formula(T)
% 
% Calculates density for specified temperature.
%
% INPUT:    T --> temperature in degree Celsius
%
% OUTPUT: rho --> density in mg/mm^3
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de



% formula and parameter from PCM_rho.m from Robert
p0 = 167.2182;
p1 = 129.9861;
p2 = 760.8218;
p3 = 0.078916;


rho = p0./(1.0+exp(p3.*(T-p1))) + p2;
rho = rho * 1e-3;  % rescale to [mg / mm^3]

return