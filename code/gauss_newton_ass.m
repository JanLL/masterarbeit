function [x_end] = gauss_newton_ass(F1_func, F2_func, x_start, lb, ub, options)
% Gauss Newton Active Set Strategy to solve inequality constrained Least
% Squares Problems
% min || F1 + J1*dx ||
% s.t.   F2 + J2*dx =  0
%        F3 + J3*dx >= 0

% Abbreviations
TOL_ineq = options.TOL_ineq;


% Some pre calculations
n = length(x_start);


idx_lb_set = find(lb ~= -inf);
idx_ub_set = find(ub ~= inf);

num_lb = length(idx_lb_set);
num_ub = length(idx_ub_set);

F3_func = @(x) [x - lb'; -(x - ub')];


F2 = F2_func(x_start);
F3 = F3_func(x_start);

m2 = length(F2);

% Check feasibility of x_start
if (sum(F2 ~= 0) > 0)
    error('x_start is not feasible w.r.t. equality constraints!')
elseif (sum(F3 < 0) > 0)
    error('x_start is not feasible w.r.t. inequality constraints!')
end

% Get initial Active Set A
A = F3 < TOL_ineq & F3 > -TOL_ineq;
A_lb = A(1:n);
A_ub = A(n+1:2*n);


x_k = x_start;
dx_norm = inf;

while (dx_norm > TOL)
    
    [F1, J1] = F1_func(x_k);
    [F2, J2] = F2_func(x_k);
    [F3, J3] = F3_func(x_k);
    
    % Build total current equality constraints from F2 and actice F3
    F_active = [F2; F3(A)];
    J_active = [J2; J3(A,:)];
    
    % Solve equality constraint LSQ subproblem
    [dx, Q1, R_bar] = gauss_newton_step_constr(F1, J1, F_active, J_active, options);
    
    % TODO: Schrittweitensteuerung "Backtracking Linesearch"
    x_kp1 = x_k + dx;
    
    % Check for violations of non-active inequality constraints
    F3_kp1 = F3_func(x_kp1);
    ineq_violation_lb = F3_kp1(1:n) < -TOL_ineq;
    ineq_violation_ub = F3_kp1(n+1:end) < -TOL_ineq;
    
    % For violations, set to value of lb, ub
    x_kp1(ineq_violation_lb) = lb(ineq_violation_lb);
    x_kp1(ineq_violation_ub) = ub(ineq_violation_ub);
    
    % Compute lagrange multipliers lambda with modified dx
    dx_mod = x_kp1 - x_k;
    lambda = R_bar \ (Q1.' * J1.' * J1 * dx_mod + Q1.' * J1.' * F1);
    
    % Remove constraints from active set with negative lambda
    lambda_aux_lb = inf*ones(n,1);
    lambda_aux_lb(A_lb == true) = lambda(m2 + (1:sum(A_lb)));
    A_lb(A_lb == true & lambda_aux_lb.' < 0) = false; % TODO: vllt mit ner Tolerance
    
    lambda_aux_ub = inf*ones(n,1);
    lambda_aux_ub(A_ub == true) = lambda(m2+sum(A_lb) + (1:sum(A_ub)));
    A_ub(A_ub == true & lambda_aux_ub.' < 0) = false; % TODO: vllt mit ner Tolerance
    
    
    % Add index of violated constraints to active set.
    A_lb(ineq_violations_lb) = true;
    A_ub(ineq_violations_ub) = true;
    
end

end

