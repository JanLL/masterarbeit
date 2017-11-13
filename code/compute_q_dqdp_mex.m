function [varargout] = compute_q_dqdp_mex(...
    p_optim_free, p_optim_estimable, p_optim_fixed, c_p_param_type, T_ref_dsc, q_dsc, ax1, ax2, ax3)
% TODO: description!

% build vector of all (free and fixed) optimization parameters
p_optim_all = zeros(1,length(p_optim_estimable));
p_optim_all(p_optim_estimable) = p_optim_free;
p_optim_all(~p_optim_estimable) = p_optim_fixed;


% c_p plot
T_plot = 30:0.01:160;
switch c_p_param_type
    case 'old_atan_formula'
        c_p_plot = c_p_formula(T_plot, p_optim_all(1:6));
    case 'fraser_suzuki'
        c_p_plot = c_p_fs(T_plot, p_optim_all);
    case 'gauss_linear_comb'
        c_p_plot = c_p_gauss_linear_comb(T_plot, p_optim_all);
end

cla(ax2); 
hold(ax2, 'on')
plot(ax2, T_plot, c_p_plot, 'DisplayName', 'c_p Simulation');
legend(ax2, 'show', 'location', 'northwest');
xlabel(ax2, 'T [degC]');
ylabel(ax2, 'c_p [mJ/(mg*K]');
drawnow;

%[residuum, Jac] = heat1D_pcm('optimization', p_optim_all); 

if (nargout == 1)
    residuum = heat1D_pcm('evaluation', p_optim_all);
    varargout{1} = residuum;
elseif (nargout == 2)    
    [residuum, Jac, T] = heat1D_pcm('optimization', p_optim_all); 
    Jac = Jac(:,p_optim_estimable);
    varargout{1} = residuum;
    varargout{2} = Jac;
end



% residuum = heat1D_pcm('optimization', p_optim_all); 

% Fuer Vergleich Fwd/Adj Sensitivities...
% [residuum, Jac_fwd, Jac_adj] = heat1D_pcm('optimization', p_optim_all);
% Jac = struct();
% Jac.fwd = Jac_fwd;
% Jac.adj = Jac_adj;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heat flux plot
q_sim = residuum + q_dsc;
cla(ax1); % q_pcm_in plot
hold(ax1, 'on')
plot(ax1, T_ref_dsc, q_sim, 'DisplayName', 'Simulation');
plot(ax1, T_ref_dsc, q_dsc, 'DisplayName', 'Measurement');
plot(ax1, T_ref_dsc, residuum, 'DisplayName', 'Residuum');
legend(ax1, 'show', 'location', 'northwest');
xlabel(ax1, 'T_{ref} [degC]');
ylabel(ax1, 'q_{pcm}^{in} [mW]');
drawnow;


% dqdp plot
if (nargout == 2)
    cla(ax3); 
    hold(ax3, 'on')
    image(ax3, Jac, 'CDataMapping', 'scaled');
    colorbar(ax3);
    title(ax3, 'dqdp')
    xlabel(ax3, 'c_p parameters');
    ylabel(ax3, 'T_{ref}');
    drawnow;
end

end

