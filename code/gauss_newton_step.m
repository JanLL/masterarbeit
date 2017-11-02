function [dx] = gauss_newton_step(F1, J1, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = size(J1,2);
m1 = size(J1,1);

if strcmp(options.decomposition,'QR')

    [Q,R] = qr(J1);

    R_bar = R(1:n,:);
    Q1 = Q(:,1:n);

    dx = (linsolve(R_bar, -Q1' * F1))';

end

if strcmp(options.decomposition,'SVD')

    [U,S,V] = svd(J1);
    
    singular_values = diag(S);
    TOL = max(size(J1))*eps(max(singular_values));
    r = sum(singular_values > TOL);  % rank computation like in matlab function rank()
        
    U1 = U(:,1:n);
    c = U1' * F1;
    
    dy = zeros(n,1);
    for i=1:r
        dy(i) = -c(i) / singular_values(i);
    end
        dy(r+1:end) = 0.;
        
    dx = (V * dy)';
        
end




end

