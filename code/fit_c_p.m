function [ output_args ] = fit_c_p()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% TODO: Bisher nur reinkopiert von calc_cp.m ... weis nicht ob ich die Fkt.
% noch brauch

% test start parameters
%fit_fct = @(x, a, b, c, d, e, f) (atan(-d*(x-a))+pi/2) .* (b*exp(-c.*(x-a).^2)) + e*x + f;
%plot(T_grid(100:end-10), fit_fct(T_grid(100:end-10), 140., 5., 0.01, 0.7, 0.01, 0.)); hold on
%plot(T_grid(100:end-10), c_p_pcm(100:end-10), 'color', 'blue');

%% fit c_p on explicit atan function
%fit_fct1 = fittype('1/(exp(d*(x-a))+1) * (b*exp(-c*(x-a)^2)) + e*x + f', ...
%    'coeff', {'a', 'b', 'c', 'd', 'e', 'f'});
% fit_fct2 = fittype('(atan(-d*(x-a))+pi/2) * (b*exp(-c*(x-a)^2)) + e*x + f', ...
%     'coeff', {'a', 'b', 'c', 'd', 'e', 'f'});
% start_param = [140., 5., 0.01, 0.7, 0.01, 0.];
% f = fit(T_grid(100:end-10), c_p_pcm(100:end-10), fit_fct2, ...
%     'Startpoint', start_param);
% 
% c_p_fit = f(T_grid(100:end-10));
% 
% %plot(T_grid(100:end-10), c_p_fit, 'color', 'red', ...
% %     'DisplayName', 'fit'); hold on
% plot(T_grid(100:end-10), c_p_pcm(100:end-10), 'color', 'blue', ...
%      'DisplayName', 'measurements');
% 
% legend(gca, 'show', 'Location', 'northwest')

%% fit c_p on BSpline

[T_ref, c_p] =  calc_cp();

index_30 = find(T_ref > 30, 1, 'first');


T_ref_knots = [find(T_ref > 30, 1, 'first'), ...
               find(T_ref > 70, 1, 'first'), ...
               find(T_ref > 100, 1, 'first'), ...
               find(T_ref > 120, 1, 'first'), ...
               find(T_ref > 130, 1, 'first'), ...
               find(T_ref > 135, 1, 'first'), ...
               find(T_ref > 140, 1, 'first'), ...
               find(T_ref > 145, 1, 'first'), ...
               find(T_ref > 150, 1, 'first'), ...
               find(T_ref > 159, 1, 'first')];

fit = csapi(T_ref(T_ref_knots), c_p(T_ref_knots));
fit_B = fn2fm(fit, 'B-');

plot(T_ref(index_30:end), c_p(index_30:end)); hold on
% plot(T_ref(index_30:end), ppval(fit, T_ref(index_30:end)));

end

