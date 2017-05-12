function [J_lin_sparse] = build_linear_matrix(N)
% [J_lin_sparse] = build_linear_matrix(N)
%  TODO AB HIER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Calculates specific heat capacity for given temperature range.
%
% INPUT:    T --> temperature in degree Celsius
%
% OUTPUT: c_p --> specific heat capacity in [mJ/(mg*K)]
%
% Author: Jan Lammel, lammel@stud.uni-heidelberg.de


J_lin_columns = ones(N, 3);

% main diagonal
J_lin_columns(2:end-1, 1) = -2.;
J_lin_columns(1, 1) = 0;
J_lin_columns(N, 1) = -1.;

% upper first diagonal
J_lin_columns(1, 2) = 0.;  % one element less than in main diagonal

% lower first diagonal
J_lin_columns(1:2, 2) = 0.;  % one element less again and one zero entry

% build actual sparse matrix for linear part
diagonals = [0, 1, -1];
J_lin_sparse = spdiags(J_lin_columns, diagonals, N, N);


end
