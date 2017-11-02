h = @(x,A,mu,var) A*exp(-(x-mu).^2/var);

A = 2.;
mu = 2.;
var = 1.;

x = (-5:0.1:5)';
y = h(x,A,mu,var);

eps = 0.05;
eta = y + (rand(length(y),1)-0.5) * 2*eps;

% plot(x, y); hold on
% plot(x, y_tilde, 'x')

A0 = 1.3;
mu0 = 1.3;
var0 = 2.;

p0 = [A0, mu0, var0];

F1_fct = @(p) h(x,p(1),p(2),p(3)) - eta;

J1_fct = @(p) [exp(-(x-p(2)).^2/p(3)), ...
                  2*p(1)*(x-p(2))/p(3) .* exp(-(x-p(2)).^2/p(3)), ...
                  p(1)*(x-p(2)).^2 / p(3)^2 .* exp(-(x-p(2)).^2/p(3))];
        
F1 = F1_fct(p0);
J1 = J1_fct(p0);     
 
% figure()
% image(J1, 'CDataMapping', 'scaled')
% colorbar

GN_options = struct;
GN_options.decomposition = 'SVD';

p = p0;
dp_norm = inf;
while dp_norm > 1e-8

    F1 = F1_fct(p);
    J1 = J1_fct(p); 
    
    dp = gauss_newton_step(F1,J1,GN_options);

    p = p + dp;
    dp_norm = norm(dp);
    fprintf('%e\n', dp_norm);
    
end


