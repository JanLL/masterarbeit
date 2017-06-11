knots = [0, 1.5, 2, 4, 4.1, 5., 6., 6.5, 7., 8., 9., 10.];
coeffs = [2., 6, -2, 1, 1, 1, 2, 1, 2];

sp = spmak(knots, coeffs);
dsp = fnder(sp);


x = linspace(0,10,2000000);

f = @(x) spval(spmak(knots, coeffs), x);

% Bspline
tic;
y = spval(sp, x);
toc

tic;
y = f(x);
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

% plot(x, spval(sp, x)); hold on;
% plot(x, spval(dsp, x)); hold on;
