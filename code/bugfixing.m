% file_src = '/home/argo/masterarbeit/simulationen-data/pcm_bspline_revMass_k-45_lsqnonlin_with_coeff_bounds/fit_data.mat';
% 
% fit_data = load(file_src);
% 
% save_path = '/home/argo/masterarbeit/fits_data/';
% 
% dsc_data_struct = fit_data.measurements.dsc_data_struct;
% index_T_dsc = fit_data.measurements.index_T_dsc;
% revMassNorm = fit_data.measurements.revMassNorm;
% p_sim = fit_data.simulation;
% optim_solverName = fit_data.optimization.solverName;
% optim_options = fit_data.optimization.optim_options;
% p_optim_start = fit_data.optimization.param_start;
% p_optim_estimable = fit_data.optimization.estimable;
% optim_con = fit_data.optimization.optim_con;
% p_optim_end = fit_data.optimization.param_end;
% optim_output = fit_data.optimization.output;
% 
% 
% save_fit(save_path, dsc_data_struct, index_T_dsc, revMassNorm, ...
%     p_sim, optim_solverName, optim_options, p_optim_start, p_optim_estimable, ...
%     optim_con, p_optim_end, optim_output);


% knots = [-10,0, 30, 50, 60, 70, 80, 90, 100, 110, 115, 122, 127, 132, 135:160, 165, 170, 200];
% coeffs = (0.05 .* [1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 5, 15, 15, 15, 15, 15, ...
%                   15, 15, 15, 15, 15, 15, 15, 15, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]);

% knots = 10:10:50;
% coeffs = 1:1:5;
%               
% n_knots = length(knots);
% n_coeffs = length(coeffs);
% 
% A_columns = ones(n_knots, 2);
% 
% % knots
% A_columns(n_knots:end,1) = 0;
% 
% % coeffs
% A_columns(1:n_knots-1,2) = -1;
% A_columns(n_knots:end,2) = 0;
% 
% diagonals = [0,1];
% A_sparse = spdiags(A_columns, diagonals, n_knots-1, n_knots+n_coeffs);
% 
% A = full(A_sparse);



string = 'blubb';

switch string
    case 'hallo'
        disp(1)
    case 'blubb'
        disp(2)
end



