function [sim_params] = get_param_sim(varargin)
%
%
% INPUT (variable):          
%             N1 --> number of spatial discretization lattice points 
%                    (Constantan)
%             N2 --> number of spatial discretization lattice points 
%                    (Crucible)
%             N3 --> number of spatial discretization lattice points (PCM)
%             L1 --> length [mm] of Constantan part.
%             L2 --> length [mm] of crucible part.
%             L3 --> length [mm] of PCM part.
%            T_0 --> initial temperature [degree celsius] everywhere.
%          T_end --> integrate until T_oven as reached T_end [degree celsius].
%      heat_rate --> rate [K/min] with which oven temperature increases.
% c_p_test_setup --> 2x1 array: specific heat capacity [mJ/(mg*K] 
%                       of [Constantan, crucible]
% rho_test_setup --> 2x1 array: density [mg/mm^3] of [Constantan, crucible]
%         lambda --> 3x1 array: heat conductivity [mW/(mm*K)] 
%                       of [Constantan, crucible, PCM]
%     c_p_sample --> 2x1 cell array: Fct Handles with unspecified
%                       parameters to eval. c_p (1) and dc_p (2) [mJ/(mg*K]
%                       of sample.
%
% OUTPUT: 
%   sim_params --> 1x2 struct with (default) simulation parameters as input 
%                  for simulate_1d. 
%                  sim_params(1): simulation with PCM
%                  sim_params(2): reference simulation without PCM

sim_params = struct();

% check for input arguments and update variables where necessary
if hasOption(varargin, 'L1'), sim_params.L1 = getOption(varargin, 'L1'); 
else sim_params.L1 = 25.; end
if hasOption(varargin, 'L2'), sim_params.L2 = getOption(varargin, 'L2'); 
else sim_params.L2 = 0.; end
if hasOption(varargin, 'L3'), sim_params.L3 = getOption(varargin, 'L3'); 
else sim_params.L3 = 1.; end

if hasOption(varargin, 'N3'), sim_params.N3 = getOption(varargin, 'N3');
else sim_params.N3 = 50; end
if sim_params.N3 <= 0, error('N3 must be greater than 0!'); end

% compute N1, N2 s.t. dx is equal everywhere as default
if hasOption(varargin, 'N1'), sim_params.N1 = getOption(varargin, 'N1'); 
else sim_params.N1 = sim_params.L1 / sim_params.L3 * sim_params.N3; end

if hasOption(varargin, 'N2'), sim_params.N2 = getOption(varargin, 'N2'); 
else sim_params.N2 = sim_params.L2 / sim_params.L3 * sim_params.N3; end

if hasOption(varargin, 'heat_rate'), sim_params.heat_rate = getOption(varargin, 'heat_rate'); 
else sim_params.heat_rate = 10; end
if hasOption(varargin, 'T_0'), sim_params.T_0 = getOption(varargin, 'T_0'); 
else sim_params.T_0 = 10.; end
if hasOption(varargin, 'T_end'), sim_params.T_end = getOption(varargin, 'T_end'); 
else sim_params.T_end = 300.; end

% src of material constants: Constantan, crucible(Al2O3) -> Wikipedia; 
% density of pcm: Robert: PCM_lambda.m
if hasOption(varargin, 'c_p_test_setup')
    sim_params.c_p_test_setup = getOption(varargin, 'c_p_test_setup'); 
else sim_params.c_p_test_setup = [0.41, 0.99]; 
end
if hasOption(varargin, 'rho_test_setup')
    sim_params.rho_test_setup = getOption(varargin, 'rho_test_setup'); 
else sim_params.rho_test_setup = [8.9, 3.94]; 
end
if hasOption(varargin, 'lambda_test_setup')
    sim_params.lambda_test_setup = getOption(varargin, 'lambda_test_setup'); 
else sim_params.lambda_test_setup = [23., 35.6, 0.96]; 
end




if hasOption(varargin, 'c_p_sample')
    c_p_sample = getOption(varargin, 'c_p_sample');
    
    if isa(c_p_sample{1}, 'function_handle') && ...
       isa(c_p_sample{2}, 'function_handle') && ...
       isa(c_p_sample{3}, 'numeric')
        sim_params.eval_c_p = c_p_sample{1};
        sim_params.eval_dc_p = c_p_sample{2};
        sim_params.c_p_params_num = c_p_sample{3};

        sim_params.get_param_c_p = ...
            @(p_optim) p_optim(1 : sim_params.c_p_params_num);
        sim_params.get_param_k = ...
            @(p_optim) p_optim(sim_params.c_p_params_num + 1);


    elseif strcmp(c_p_sample{1}, 'B-') && ...
           isa(c_p_sample{2}, 'numeric')
        sim_params.c_p_params_num = c_p_sample{2};
        c_p_params_num = sim_params.c_p_params_num;

        sim_params.get_param_c_p_knots = ...
            @(p_optim) p_optim(1 : c_p_params_num(1));
        sim_params.get_param_c_p_coeffs = ...
            @(p_optim) p_optim(c_p_params_num(1)+1 : sum(c_p_params_num(1:2)));
        sim_params.get_param_c_p = ...
            @(p_optim) cat(2, sim_params.get_param_c_p_knots(p_optim), ...
                              sim_params.get_param_c_p_coeffs(p_optim));
        sim_params.get_param_k = ...
            @(p_optim) p_optim(sum(c_p_params_num(1:2))+1);

        sim_params.eval_c_p = ...
            @(T, p) spval(spmak(sim_params.get_param_c_p_knots(p), ...
                                sim_params.get_param_c_p_coeffs(p)), T);
        sim_params.eval_dc_p = ...
            @(T, p) spval(fnder(spmak(sim_params.get_param_c_p_knots(p), ...
                                sim_params.get_param_c_p_coeffs(p))), T);
    else
        error(['Option c_p_sample must be either ' ...
              '{@cp_fct, @dcp_fct, [num_params]} for function handles' ...
              'or {''B-'', [#knots, #coeffs]} for BSplines.']);
    end
    
else
    syms T;
    p = sym('p', [6, 1]);
    dc_p = matlabFunction(diff(c_p_formula(T, p), T), 'Vars', [T;p]);
    dc_p = @(T,p) dc_p(T,p(1),p(2),p(3),p(4),p(5),p(6));
    
    sim_params.eval_c_p = @c_p_formula;
    sim_params.eval_dc_p = dc_p;
    
    sim_params.c_p_params_num = 6;
    sim_params.get_param_c_p = ...
        @(p_optim) p_optim(1 : sim_params.c_p_params_num);
    sim_params.get_param_k = ...
        @(p_optim) p_optim(sim_params.c_p_params_num + 1);
end

    
    
sim_params(2) = deal(sim_params(1));
sim_params(2).N3 = 0;
sim_params(2).L3 = 0.;




end


