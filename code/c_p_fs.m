function [c_p] = c_p_fs(T, p)
% [c_p] = c_p_fs(T)
% 
% Calculates specific heat capacity for given temperature range using model
% of a Fraser-Suzuki Peak.
%
% INPUT:            T --> temperature in [degree Celsius]
%          c_p_params --> parameters of function, 6x1 array
%
% OUTPUT: c_p --> specific heat capacity in [mJ/(mg*K)]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de

h = p(1);
r = p(2);
wr = p(3);
sr = p(4);
z = p(5);
m = p(6);
b = p(7);


% determine the nonzero indices
nonzeroIDX = ( T < z - (wr.*sr)./(sr^2-1) );
Tnonzero = T(nonzeroIDX);

% first output argument: nominal values
c_p = m*T + b;
c_p(nonzeroIDX) = c_p(nonzeroIDX) + h*exp(-log(r)/log(sr)^2 * (log(1+(Tnonzero-z)*(sr^2-1)/(wr*sr))).^2);


return