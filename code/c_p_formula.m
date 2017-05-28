function [c_p] = c_p_formula(T, p)
% [c_p] = c_p_formula(T)
% 
% Calculates specific heat capacity for given temperature range.
%
% INPUT:            T --> temperature in [degree Celsius]
%          c_p_params --> parameters of function, 6x1 array
%
% OUTPUT: c_p --> specific heat capacity in [mJ/(mg*K)]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de



    

c_p = (atan(-p(4)*(T-p(1)))+pi/2) .* (p(2)*exp(-p(3)*(T-p(1)).^2)) + p(5)*T + p(6);
% atan works as activation function to get an asymmetric peak

return