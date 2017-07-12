N1 = 4;
N3 = 6;

N = N1+N3;


J_lin = build_linear_matrix(N);

J_lin(N1,N1-1) = 66;
J_lin(N1,N1) = 5;
J_lin(N1,N1+1) = 77;

full(J_lin)
