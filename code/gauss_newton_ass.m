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

% Check feasibility of x_start
if (sum(F2 ~= 0) > 0)
    error('x_start is not feasible w.r.t. equality constraints!')
elseif (sum(F3 < 0) > 0)
    error('x_start is not feasible w.r.t. inequality constraints!')
end

% Get initial Active Set A
A = F3 < TOL_ineq & F3 > -TOL_ineq;


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
    [dx, lambda] = gauss_newton_step_constr(F1, J1, F_active, J_active, options);
    
    % TODO: Schrittweitensteuerung "Backtracking Linesearch"
    x_kp1 = x_k + dx;
    
    % Check for violations of non-active inequality constraints
    F3_kp1 = F3_func(x_kp1);
    ineq_violation_lb = F3_kp1(1:n) < -TOL_ineq;
    ineq_violation_ub = F3_kp1(n+1:end) < -TOL_ineq;
    
    % For violations, set to value of lb, ub
    x_kp1(ineq_violation_lb) = lb(ineq_violation_lb);
    x_kp1(ineq_violation_ub) = ub(ineq_violation_ub);
    
    % Add corresponding index to active set.
    A_lb(ineq_violations_lb) = true;
    A_ub(ineq_violations_ub) = true;
    
    % Recompute lambda with modified dx
    % TODO!
    
    % Remove constraints from active set with negative lambda
    % TODO!
    
end

end

