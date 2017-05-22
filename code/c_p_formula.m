function [c_p] = c_p_formula(T, heat_rate)
% [c_p] = c_p_formula(T)
% 
% Calculates specific heat capacity for given temperature range.
%
% INPUT:            T --> temperature in [degree Celsius]
%           heat_rate --> Shape (blurring) of c_p curve depends on the
%                         the heat rate (0 for approx delta-peak).
%
% OUTPUT: c_p --> specific heat capacity in [mJ/(mg*K)]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de


% values from fit of measurement for beta=10K/min
p0 = 144.0009;
p1 = 4.1036;
p2 = 0.0039;
p3 = 1.4217; p3=0;
p4 = 0.0078;
p5 = 1.5325;

if heat_rate == 0
    p1 = p1 * 5.;
    p2 = p2 + 0.1;
end
    

c_p = (atan(-p3*(T-p0))+pi/2) .* (p1*exp(-p2*(T-p0).^2)) + p4*T + p5;
% atan works as activation function to get an asymmetric peak

return