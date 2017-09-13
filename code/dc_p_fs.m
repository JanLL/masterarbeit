function [dc_p] = dc_p_fs(T, p)
% [dc_p] = dc_p_fs(T)
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
b = p(6);

% determine the nonzero indices
nonzeroIDX = ( T < z - (wr.*sr)./(sr^2-1) );
Tnonzero = T(nonzeroIDX);

log_arg = 1+(Tnonzero-z)*(sr^2-1)/(wr*sr);
exp_arg = -log(r)/(log(sr)^2) * (log(log_arg)).^2;


dc_p = zeros(size(T));


dc_p(nonzeroIDX) = -2*log(r)/log(sr)^2 * (sr^2-1)/(wr*sr) * log(log_arg) ./ log_arg .* h.*exp(exp_arg);


return