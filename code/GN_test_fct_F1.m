function [F1, J1] = GN_test_fct_F1(p, x, eta)
%

h = @(x,p) p(1)*exp(-(x-p(3)).^2/p(2));

F1_eval_fct = @(p) h(x,p) - eta;

J1_eval_fct = @(p) [exp(-(x-p(3)).^2 / p(2)), ...
               p(1) * (x-p(3)).^2 / p(2)^2 .* exp(-(x-p(3)).^2 / p(2)), ...
               2*p(1) * (x-p(3))/p(2) .* exp(-(x-p(3)).^2 / p(2))];
                
F1 = F1_eval_fct(p);
J1 = J1_eval_fct(p);


% figure(1); clf; hold on
% plot(x, y);
% plot(x, eta, 'x');
% plot(x, F1);

return

end

