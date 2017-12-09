function [theta, C] = compute_confidence_interval(F1, J1, alpha, scaling)
% Currently just implemented for unconstrained case!
% 
% INPUT:
%   F1    -> Residuum column vector
%   J1    -> Jacobian of F1 w.r.t. opt. parameters
%   alpha -> Level of significance, e.g. 0.05
% scaling -> Vector of scaling vectors.
%
% OUTPUT:
%   theta -> Cuboid confidence interval
%       C -> Covariance matrix


n_mp = size(J1,1);
n_p = size(J1,2);


% Compute Covariance matrix
[~,Sigma,V] = svd(full(J1));
    
Sigma_inv_square = diag(diag(Sigma).^(-2));

C = V * Sigma_inv_square * V.';


% Get measurement error factor b^2
b = sqrt(F1.' * F1 / (n_mp - n_p));


% Compute (1-alpha) quantile of Fisher-Distr.
gamma = n_p * finv(1-alpha, n_p, n_mp - n_p);


% Compute final cuboid confidence interval theta
theta = b * sqrt(gamma * diag(C).');

end

