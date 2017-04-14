function [c_p] = c_p_formula(T)

% values from fit for beta=10K/min
p0 = 144.0394;
p1 = 100.6760;
p2 = 0.0039;
p3 = 1.1780;
p4 = 0.0730;
p5 = 11.8995;

c_p = (atan(-p3*(T-p0))+pi/2) .* (p1*exp(-p2*(T-p0).^2)) + p4*T + p5;

return