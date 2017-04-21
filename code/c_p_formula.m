function [c_p] = c_p_formula(T)

% values from fit for beta=10K/min
p0 = 144.0009;
p1 = 4.1036;
p2 = 0.0039;
p3 = 1.4217;
p4 = 0.0078;
p5 = 1.5325;

c_p = (atan(-p3*(T-p0))+pi/2) .* (p1*exp(-p2*(T-p0).^2)) + p4*T + p5;

return