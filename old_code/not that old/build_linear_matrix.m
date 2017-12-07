function [J_lin_sparse] = build_linear_matrix(N)
% [J_lin_sparse] = build_linear_matrix(N)
% 
% Builds the sparse matrix of the linear part of the heat equation:
%
% |0   0  0 .... 0|
% |1  -2  1 0 .. 0|
% |0   1 -2 1 0..0|
% |...............|
% |0...0  1  -2  1|
% |0...0  0   1 -1|
%
% INPUT:             N --> number of spatial lattice points.
%
% OUTPUT: J_lin_sparse --> sparse matrix of the linear part of the 
%                          heat equation (see above structure)
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


return
