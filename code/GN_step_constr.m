function [dx, Q1, R_bar] = GN_step_constr(F1, J1, F2, J2, options)
% 


n = size(J1, 2);
m1 = size(J1, 1);
m2 = size(J2, 1);

if (length(F1) ~= m1 || size(J1,2) ~= n || size(J2,2) ~= n || length(F2) ~= m2)
    error('Dimensions of input do not fit!\n')
end
    
[Q, R] = qr(J2.');

R_bar = R(1:m2,:);
Q1 = Q(:,1:m2);

dy1 = - R_bar.' \ F2;

A = J1 * Q;
A1 = A(:,1:m2);
A2 = A(:,m2+1:end);

F1_new = F1 + A1*dy1;
J1_new = A2;

dy2 = GN_step(F1_new, J1_new, options);

dy = [dy1; dy2];


dx = Q * dy;


end

