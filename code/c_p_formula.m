function [c_p] = c_p_formula(T)

% values from fit for beta=10K/min
p0 = 140.;
p1 = 100.;
p2 = 0.01;
p3 = 0.7;
p4 = 0.01;
p5 = 15.;

c_p = (atan(-p3*(T-p0))+pi/2) .* (p1*exp(-p2*(T-p0).^2)) + p4*T + p5;

return