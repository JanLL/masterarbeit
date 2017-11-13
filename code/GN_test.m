addpath('/home/argo/')


J1=[[1,1,1,1];[1,3,1,1];[1,-1,3,1];[1,1,1,3];[1,1,1,-1]];
F1=-[2,1,6,3,1]';
J2=[[1,1,1,-1];[1,-1,1,1];[1,1,-1,1]];
F2=-[2,3,-1]';


fprintf('Loesung von Uli\n');
[dx_Uli,lambda_Uli,res1_Uli] = clls(J1,J2,F1,F2);
[dy_Uli,res2_Uli] = lls(J1,F1);

% dy_Uli'
% res2_Uli'

dx_Uli'
lambda_Uli'
% res1_Uli'


fprintf('Meine Loesung\n');
GN_options = struct;
GN_options.decomposition = 'SVD';

[dx] = GN_step(F1, J1, GN_options);


[dx, lambda] = GN_step_constr(F1, J1, F2, J2, GN_options);
dx'
lambda'





