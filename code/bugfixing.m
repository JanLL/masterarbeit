knots = [-50, 60., 90., 120., 130, 135, 140, 145, 150, 160, 220];
coeffs = [1, 1, 1, 20, 1, 1, 1];


sp = spmak(knots, coeffs);
dsp = fnder(sp);

x = 30:0.05:160;
 
% Bspline
tic;
y = spval(sp, x);
toc

% % Normales Polynom (quadratisch)
% coeffs = [2,1,1];
% tic;
% polyval(coeffs, x);
% toc
% 
% % c_p Auswertung bisher
% p_optim_start = [115., ...
%                  25., ...
%                  0.01, ...
%                  -1.42, ...
%                  0.018, ...
%                  5.4, ...
%                  0.3];
% 
% c_p_formula_expl = @(T) c_p_formula(T, p_optim_start);       
% tic;
% c_p_formula_expl(x);
% toc

%dsp = fnder(sp);

plot(x, spval(sp, x)); hold on;
plot(x, spval(dsp, x)); hold on;
