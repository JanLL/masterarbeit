%% Test with lb, ub inequality constraints using Active Set Strategy
% Synthesize measurement data
A = 2.;
var = 1.;
mu = 2.;
p_true = [A, var, mu];  % true variables

h = @(x,p) p(1)*exp(-(x-p(3)).^2/p(2));

x = (-5:0.01:5).';
y = h(x, p_true);

eps = 0.05;
eta = y + (rand(length(y),1)-0.5) * 2*eps;

% Optimization Start values
A0 = 1.25;
var0 = 1.5;
mu0 = 1.5;

p0 = [A0; var0; mu0];


GN_options = struct;
GN_options.decomposition = 'SVD';
GN_options.TOL_ineq = 1e-8;
GN_options.TOL_dx_norm = 1e-200;
GN_options.TOL_t_k = 1e-200;
GN_options.max_iterations = 1000;
GN_options.t_k_start = 0.3;

lb = [-inf, -inf, -inf];
ub = [inf, inf, inf];

F1_func = @(p) GN_test_fct_F1(p, x, eta);
F2_func = @(p) GN_test_fct_F2(p);


% Test Active Set Strategy
p_end = GN_ass(F1_func, F2_func, p0, lb, ub, GN_options);

p_end'


%% Working example without any constraints

% Synthesize measurement data
A = 2.;
var = 1.;
mu = 2.;
p_true = [A, var, mu];  % true variables

h = @(x,p) p(1)*exp(-(x-p(3)).^2/p(2));

x = (-5:0.1:5).';
y = h(x, p_true);

eps = 0.05;
eta = y + (rand(length(y),1)-0.5) * 2*eps;

% Optimization Start values
A0 = 2.1;
var0 = 1.1;
mu0 = 2.1;

p0 = [A0; var0; mu0];

p = p0;

GN_options = struct;
GN_options.decomposition = 'SVD';

dp_norm = inf;
while dp_norm > 1e-8

    [F1, J1] = GN_test_fct_F1(p, x, eta);
    
    
    dp = GN_step(F1, J1, GN_options);

    p = p + dp;
    dp_norm = norm(dp);
    fprintf('%e\n', dp_norm);
    
end

p'


