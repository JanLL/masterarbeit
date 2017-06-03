function [ output_args ] = fit_c_p()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% TODO: Bisher nur reinkopiert von calc_cp.m ... weis nicht ob ich die Fkt.
% noch brauch

% test start parameters
%fit_fct = @(x, a, b, c, d, e, f) (atan(-d*(x-a))+pi/2) .* (b*exp(-c.*(x-a).^2)) + e*x + f;
%plot(T_grid(100:end-10), fit_fct(T_grid(100:end-10), 140., 5., 0.01, 0.7, 0.01, 0.)); hold on
%plot(T_grid(100:end-10), c_p_pcm(100:end-10), 'color', 'blue');

%% fit c_p on function
%fit_fct1 = fittype('1/(exp(d*(x-a))+1) * (b*exp(-c*(x-a)^2)) + e*x + f', ...
%    'coeff', {'a', 'b', 'c', 'd', 'e', 'f'});
fit_fct2 = fittype('(atan(-d*(x-a))+pi/2) * (b*exp(-c*(x-a)^2)) + e*x + f', ...
    'coeff', {'a', 'b', 'c', 'd', 'e', 'f'});
start_param = [140., 5., 0.01, 0.7, 0.01, 0.];
f = fit(T_grid(100:end-10), c_p_pcm(100:end-10), fit_fct2, ...
    'Startpoint', start_param);

c_p_fit = f(T_grid(100:end-10));

%plot(T_grid(100:end-10), c_p_fit, 'color', 'red', ...
%     'DisplayName', 'fit'); hold on
plot(T_grid(100:end-10), c_p_pcm(100:end-10), 'color', 'blue', ...
     'DisplayName', 'measurements');

legend(gca, 'show', 'Location', 'northwest')


coeffvalues(f)

end

