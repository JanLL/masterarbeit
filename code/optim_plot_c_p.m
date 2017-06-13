function stop = optim_plot_c_p(x, optimValues, state, eval_c_p, c_p_meas, ...
                               get_param_c_p, p_optim_estimable, ...
                               p_optim_fixed, ax2)
% Plots optimized and measured c_p curves.
%
% INPUTS:
%           x --> 1D array of optimization variables.
% optimValues --> Struct containing data from the current iteration.
%       state --> current state of the algorithm: ['init', 'interrupt', ...
%                 'iter', 'done'].
%   
%    eval_c_p --> FctHandle to evaluate c_p from optimization.
%    c_p_meas --> Array of measurements of c_p. (:,1) contains temp. and 
%                 (:,2) contains c_p values.
% get_param_c_p --> FctHandle to get parameters of c_p function given all
%                   (including fixed) optimization variables.
% p_optim_estimable --> logical vector of free/fixed opt. variables.
% p_optim_fixed --> vector of opt. variable values which are fixed.
%
% OUTPUTS:
%       stop --> logical value whether to stop alg. or not.

stop = false;

cla(ax2);

p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = x;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

T = c_p_meas(1,1):0.05:c_p_meas(end,1);
plot(ax2, T, eval_c_p(T, get_param_c_p(p_optim_all)), 'DisplayName', 'Optimization'); hold on
plot(ax2, c_p_meas(:,1), c_p_meas(:,2), 'DisplayName', 'Measurement');

legend(ax2, 'show', 'location', 'northwest');
xlabel(ax2, 'T_{ref}');
ylabel(ax2, 'c_p');

drawnow;

end

