function [dc_p] = dc_p_formula(T)


% values from fit for beta=10K/min
p0 = 140.;
p1 = 100.;
p2 = 0.01;
p3 = 0.7;
p4 = 0.01;
p5 = 15.;

% previously computed via symbolic diff()  
dc_p = p4 - (p1*p3*exp(-p2*(T - p0).^2)) ./ (p3^2*(T - p0).^2 + 1) - ...
          p1*p2*exp(-p2*(T - p0).^2) .* (pi/2 - atan(p3*(T - p0))) .* ...
          (2*T - 2*p0);

return