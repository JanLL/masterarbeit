function [y] = horner_eval(x, coeffs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n = length(coeffs);
m = length(x);

y = coeffs(1) * ones(1,m);
for i = 2:n
    y = y .* x + coeffs(i);
end




end

