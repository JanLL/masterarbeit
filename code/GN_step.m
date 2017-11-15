function [dx] = GN_step(F1, J1, options)
% Compute one Gauss-Newton step in the linearized LSQ problem
% min || F1 + J1*dx ||_2^2
%
% INPUT:
%   F1      -> residuum vector to minimize
%   J1      -> Jacobian dF1/dx
%   options -> options struct to choose decomposition type
%
% OUTPUT:
%   dx      -> Gauss-Newton step of optimization variables.

n = size(J1,2);
m1 = size(J1,1);

if strcmp(options.decomposition,'QR')

    [Q,R] = qr(J1);

    R_bar = R(1:n,:);
    Q1 = Q(:,1:n);

    dx = linsolve(R_bar, -Q1' * F1);

end

if strcmp(options.decomposition,'SVD')

    [U,S,V] = svd(J1);
    
    singular_values = diag(S);
    
    % rank r computation like in matlab function rank()
    TOL = max(size(J1))*eps(max(singular_values));
    r = sum(singular_values > TOL);  
        
    U1 = U(:,1:n);
    c = U1' * F1;
    
    dy = zeros(n,1);
    for i=1:r
        dy(i) = -c(i) / singular_values(i);
    end
        dy(r+1:end) = 0.;
        
    dx = V * dy;
        
end




end

