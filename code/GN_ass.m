function [x_end, optim_output] = GN_ass(F1_func, F2_func, x_start, lb, ub, options)
% Gauss Newton Active Set Strategy to solve inequality constrained Least
% Squares Problems
% min || F1 + J1*dx ||
% s.t.   F2 + J2*dx =  0
%        F3 + J3*dx >= 0
%
% INPUT: 
%   F1_func  ->  Function handle to evaluate objective function for
%                residuum and Jacobian.
%   F2_func  ->  Function handle to evaluate equality constraints function
%                for value and Jacobian. 
%   x_start  ->  Start values of optimization variables.
%   lb       ->  Lower bound of optimization variables.
%   ub       ->  Upper bound of optimization variables.
%   options  ->  Options struct
%
% OUTPUT:
%   x_end    ->  End values of optimization variables after optimization.

tic;  % start optimization duration time measurement

% Abbreviations
TOL_ineq = options.TOL_ineq;
TOL_dx_norm = options.TOL_dx_norm;
TOL_NOC1 = options.TOL_NOC1;

TOL_t_k = options.TOL_t_k;
max_iterations = options.max_iterations;

F3_func = @(p) GN_eval_ineq_constraints(p, lb, ub);

F1_norm_vec = [];
dx_norm_vec = [];
t_k_vec     = [];
NOC1_vec    = [];


% Some pre calculations
F2 = F2_func(x_start);
F3 = F3_func(x_start);


n = length(x_start);
m2 = length(F2);
m3 = 2*n;  % fixed for just lb, ub of opt variables x

% Check feasibility of x_start
if (sum(F2 ~= 0) > 0)
    error('x_start is not feasible w.r.t. equality constraints!')
elseif (sum(F3 < 0) > 0)
    error('x_start is not feasible w.r.t. inequality constraints!')
end

% Get initial Active Set A
active_lb = F3(1:n) < TOL_ineq;
active_ub = F3(n+1:2*n) < TOL_ineq;

A = [active_lb; active_ub].';
A_old = A;

x_k = x_start;
dx_norm = inf;
NOC1 = inf;

t_k = options.t_k_start;  % initial stepsize

i = 1;
while (NOC1 > TOL_NOC1 && dx_norm > TOL_dx_norm && t_k > TOL_t_k && i <= max_iterations)
        
    t_k = 1.;
    
    if (i == 1)
        [F1, J1] = F1_func(x_k);
        [F2, J2] = F2_func(x_k);
        [F3, J3] = F3_func(x_k);

        % Build total current equality constraints from F2 and actice F3
        F2_active = [F2; F3(A)];
        J2_active = [J2; J3(A,:)];

        % Solve equality constraint LSQ subproblem
        [dx, Q1, R_bar] = GN_step_constr(F1, J1, F2_active, J2_active, options);    
        lambda = R_bar \ (Q1.' * (J1.' * J1) * dx + Q1.' * J1.' * F1);
        
        F1_norm = norm(F1);
        if (isempty(lambda) || isempty(J2_active))
            NOC1 = norm(2*J1.'*F1);
        else
            NOC1 = norm(2*J1.'*F1 - J2_active.'*lambda);
        end
        
        
        F1_norm_vec = [F1_norm_vec; F1_norm];
        NOC1_vec = [NOC1_vec; NOC1];

    end
        
    % Step size control backtracking linesearch with Armojo strategy
    c = 0.;  % parameter of Armijo Strategy (0, 0.5) for steepness
    d = 0.8;  % parameter of t_k decrease rate
    
    % some pre-computations
    F1_norm_k = norm(F1);
    J1J1 = J1.' * J1;
    dxdx = dx.' * dx;
    
        
    % cut dx such that lower/upper bounds are not violated
    A_tmp = A | A_old;
    dx_cut_factor = max(t_k * [dx(~A_tmp(1:n)); dx(~A_tmp(n+1:2*n))] ./ ...
        ([lb(~A_tmp(1:n)).'; ub(~A_tmp(n+1:2*n)).'] - [x_k(~A_tmp(1:n)); x_k(~A_tmp(n+1:2*n))]));
    
    if (dx_cut_factor > 1)
        fprintf('Cut off: %1.2f\n', dx_cut_factor);
        dx_tmp = dx / dx_cut_factor;
    else
        dx_tmp = dx;
    end
    
    x_kp1 = x_k + t_k * dx_tmp;
    
    % Line Search with try-block to check for integration errors
    try
        F1_kp1 = F1_func(x_kp1);
        F1_norm_kp1 = norm(F1_kp1);
    catch
        F1_norm_kp1 = inf;
    end
        
    if F1_norm_kp1 < F1_norm_k - c*t_k * (dxdx + dx.' * J1J1 * dx)
        t_k = min(t_k / d, 1.);
    else
        while (F1_norm_kp1 >= F1_norm_k - c*t_k * (dx.'*dx) && t_k > TOL_t_k)
            
            x_kp1 = x_k + t_k * dx_tmp;
            try  % Catch error if integration not successful (e.g. when c_p negative)
                F1_kp1 = F1_func(x_kp1);
            catch ME
                fprintf('Error occured at integration: %s\n', ME.identifier);
                t_k = t_k * d;
                
                if t_k > TOL_t_k
                    continue
                else
                    x_kp1 = x_k;
                    break
                end
            end    
            F1_norm_kp1 = norm(F1_kp1);        
            t_k = t_k * d;    
            
            if (t_k < TOL_t_k)
                x_kp1 = x_k;
                break
            end
            
        end
    end
    
    
    x_k = x_kp1;
    
    [F1, J1] = F1_func(x_k);
    [F2, J2] = F2_func(x_k);
    [F3, J3] = F3_func(x_k);

    % Update active set
    active_lb = F3(1:n) < 1e-8;
    active_ub = F3(n+1:2*n) < 1e-8;
    
    A_old = A;
    A = [active_lb; active_ub].';
        
    % Build total current equality constraints from F2 and active F3
    F_active = [F2; F3(A)];
    J2_active = [J2; J3(A,:)];

    % Solve equality constraint LSQ subproblem
    if (NOC1_vec(end) > 0.1 || true)
        [dx, Q1, R_bar] = GN_step_constr(F1, J1, F_active, J2_active, options);    
        lambda = R_bar \ (Q1.' * (J1.' * J1) * dx + Q1.' * J1.' * F1);
    else
        % Gradient direction for testing (currently set inactive!)
        dx = - 2 * J1.' * F1;
        lambda = [];
    end
    
    % Save one constraint from active set with negative lambda to remove
    % later
    lambda_aux = inf*(ones(2*n,1));
    lambda_aux(reshape(A,1,[]) == true) = lambda(m2 + (1:sum(sum(A))));
    
    idx_lambda_aux_neg = find(lambda_aux < 0);
    
    if (isempty(idx_lambda_aux_neg))
        idx_remove_temp = [];
    else
        % Choose one negative lambda randomly to remove later
        idx_remove_temp = randi(length(idx_lambda_aux_neg));
    end
    idx_remove_constraint = idx_lambda_aux_neg(idx_remove_temp);    
    
    % Remove constraint with negative lambda from previous calculation
    A(idx_remove_constraint) = false;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dx_norm = sqrt(dx.' * dx);
    F1_norm = norm(F1);
    if (isempty(lambda) || isempty(J2_active))
        NOC1 = norm(2*J1.'*F1);
    else
        NOC1 = norm(2*J1.'*F1 - J2_active.'*lambda);
    end
    
    % Save optimization process variables
    F1_norm_vec = [F1_norm_vec; F1_norm];
    dx_norm_vec = [dx_norm_vec; dx_norm];
    t_k_vec     = [t_k_vec; t_k];
    NOC1_vec    = [NOC1_vec; NOC1];
    
    fprintf('\nIteration: %d\tF1_norm: %1.3e\tdx_norm: %1.3e\tt_k = %1.3e\tNOC1: %1.3e\n', i, F1_norm, dx_norm, t_k, NOC1);
    
    for j=0:ceil(length(x_k)/10)-1
        fprintf('param number:\t');
        fprintf(repmat('%d\t',1,10),(1:10)+j*10 );
        fprintf('\n');

        fprintf('param value:\t');
        index = (1:10)+j*10;
        index = index(index <= length(x_k));
        fprintf(repmat('%1.3g\t',1,10), x_k(index));
        fprintf('\n');
        
        % TODO: Active Set!
    end
    
    i = i+1;
    
end


x_end = x_k;

optim_duration = toc;  % end optimization duration time measurement


optim_output = struct();
optim_output.progress_F1_norm = F1_norm_vec;
optim_output.progress_dx_norm = dx_norm_vec;
optim_output.progress_t_k = t_k_vec;
optim_output.progress_NOC1 = NOC1_vec;
optim_output.optim_duration = optim_duration;


end

