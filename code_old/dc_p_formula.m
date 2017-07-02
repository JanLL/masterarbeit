function [dc_p] = dc_p_formula(T, p)
% [dc_p] = dc_p_formula(T)
% 
% Calculates d(specific heat capacity)/dT for given temperature range.
%
% INPUT:         T --> temperature in [degree Celsius]
%       c_p_params --> parameters of function, 6x1 array
%
% OUTPUT: dc_p --> d(specific heat capacity)/dT in [mJ/(mg*K^2)]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de


% previously computed derivative from c_p_formula.m via symbolic diff()  
dc_p = p(5) - (p(2)*p(4)*exp(-p(3)*(T - p(1)).^2)) ./ (p(4)^2*(T - p(1)).^2 + 1) - ...
          p(2)*p(3)*exp(-p(3)*(T - p(1)).^2) .* (pi/2 - atan(p(4)*(T - p(1)))) .* ...
          (2*T - 2*p(1));

return