function [dc_p] = dc_p_formula(T)
% [dc_p] = dc_p_formula(T)
% 
% Calculates d(specific heat capacity)/dT for given temperature range.
%
% INPUT:     T --> temperature in [degree Celsius]
%
% OUTPUT: dc_p --> d(specific heat capacity)/dT in [mJ/(mg*K^2)]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de

% values from fit for beta=10K/min
p0 = 144.0009;
p1 = 4.1036;
p2 = 0.0039;
p3 = 1.4217;
p4 = 0.0078;
p5 = 1.5325;

% previously computed derivative from c_p_formula.m via symbolic diff()  
dc_p = p4 - (p1*p3*exp(-p2*(T - p0).^2)) ./ (p3^2*(T - p0).^2 + 1) - ...
          p1*p2*exp(-p2*(T - p0).^2) .* (pi/2 - atan(p3*(T - p0))) .* ...
          (2*T - 2*p0);

return