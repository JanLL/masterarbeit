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

scale_h  = 14.;
scale_r  = 2.;
scale_wr = 10.7;
scale_sr = 0.705;
scale_z  = 129.;
scale_m  = 0.00789;
scale_b  = 1.69;

h = p(1) * scale_h;
r = p(2) * scale_r;
wr = p(3) * scale_wr;
sr = p(4) * scale_sr;
z = p(5) * scale_z;
m = p(6) * scale_m;
b = p(7) * scale_b;


% determine the nonzero indices
nonzeroIDX = ( T < z - (wr.*sr)./(sr^2-1) );
Tnonzero = T(nonzeroIDX);

% first output argument: nominal values
c_p = 0.01*m*T + b;
c_p(nonzeroIDX) = c_p(nonzeroIDX) + h*exp(-log(r)/log(sr)^2 * (log(1+(Tnonzero-z)*(sr^2-1)/(wr*sr))).^2);


return