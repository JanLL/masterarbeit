% testing of rational polynome function

P_coeffs = [0., 0.1, 0, 20.,  0, 1, 15, 3];
Q_coeffs = [1., 0.05, 0, -20, 1, -10, 15, 3];


x = 0.1:0.01:20;

y = polyval(P_coeffs, x) ./ polyval(Q_coeffs, x);

figure(1)
plot(x, y)



return

% testing of old own c_p parametrization with arctan

syms T;
p = sym('p', [6, 1]);
dc_p_formula = matlabFunction(diff(c_p_formula(T, p), T), 'Vars', [T;p]);
dc_p_formula = @(T,p) dc_p_formula(T,p(1),p(2),p(3),p(4),p(5),p(6));

dc_p_derivation = @(T,p) -(p(2)*p(4)*exp(-p(3)*(T-p(1)).^2))./(1+(p(4)*(T-p(1))).^2) ...
                         - 2*p(2)*p(3)*(T-p(1)).*exp(-p(3)*(T-p(1)).^2) ...
                         .* (atan(-p(4)*(T-p(1)))+pi/2) + p(5);


p = [130., 10, 0.01, 0.5, 0.003, 2];

T = 30:0.01:160;
c_p = c_p_formula(T, p);
dc_p = dc_p_formula(T, p);

dc_p_der = dc_p_derivation(T, p);

figure(1)
clf;
plot(T, c_p); hold on
plot(T, dc_p)
plot(T, dc_p_der, '--')

return











% testing of frazer-suzuki parametrization

r  =   20;
h  =  40;
z  = 120.0;
wr =  25;
sr =   0.4;
x  = 3:0.2:250;


[y, dy] = frazersuzuki(r, h, z, wr, sr, x);


clf;
plot(x, y); hold on
plot(x, dy)
