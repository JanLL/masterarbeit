function dae_test()
        
startup;
        
num_T = N;
num_p = N-1;

% DAE system of the form M*y' = f(y,t)
M_aux = ones(num_T+num_p,1);
M_aux(end-num_p+1:end) = 0;
M = diag(M_aux);




t_domain = linspace(0, 10);
y0 = [0.];

%ode15s(f, t_domain, y0);

