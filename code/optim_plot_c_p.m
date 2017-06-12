function stop = optim_plot_c_p(x, optimValues, state, eval_c_p, c_p_meas, ...
                               get_param_c_p, p_optim_estimable, p_optim_fixed)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

stop = false;

fig2 = figure(2);
clf(fig2);

p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = x;
p_optim_all(~p_optim_estimable) = p_optim_fixed;

T = c_p_meas(1,1):0.05:c_p_meas(end,1);
plot(T, eval_c_p(T, get_param_c_p(p_optim_all)), 'DisplayName', 'Optimization'); hold on
plot(c_p_meas(:,1), c_p_meas(:,2), 'DisplayName', 'Measurement');

legend('show', 'location', 'northwest');
xlabel('T_{ref}');
ylabel('c_p');

drawnow;

end

