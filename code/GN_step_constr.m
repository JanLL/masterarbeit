function [dx, Q1, R_bar] = GN_step_constr(F1, J1, F2, J2, options)
% Compute one Gauss-Newton step in the linearized eq-constrained LSQ problem
% min  || F1 + J1*dx ||_2^2
% s.t.    F2 + J2*dx = 0
%
% INPUT:
%   F1      -> residuum vector to minimize
%   J1      -> Jacobian dF1/dx
%   F2      -> eq constraints vector
%   J2      -> Jacobian dF2/dx
%   options -> options struct
%
% OUTPUT:
%   dx      -> Gauss-Newton step of optimization variables.
%   Q1      -> Q-part of QR decomposition of J2.
%   R_bar   -> R-part of QR decomposition of J2.



n = size(J1, 2);
m1 = size(J1, 1);
m2 = size(J2, 1);

% Check if input dimensions are correct
if (length(F1) ~= m1 || size(J1,2) ~= n || size(J2,2) ~= n || length(F2) ~= m2)
    error('Dimensions of input do not fit!\n')
end
    
% QR decomposition of J2^T
[Q, R] = qr(J2.');

R_bar = R(1:m2,:);
Q1 = Q(:,1:m2);

% Solution of constraint
dy1 = - R_bar.' \ F2;

% Build with solution dy1 new unconstrained subproblem
A = J1 * Q;
A1 = A(:,1:m2);
A2 = A(:,m2+1:end);

F1_new = F1 + A1*dy1;
J1_new = A2;

% Solve new unconstrained subproblem with existing code
dy2 = GN_step(F1_new, J1_new, options);

dy = [dy1; dy2];

dx = Q * dy; % final solution


end

