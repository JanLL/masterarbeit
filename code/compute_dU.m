function [dU_output] = compute_dU(varargin)
% TODO: Update description!
%
% INPUT:
%    fit_data_input --> either struct with saved fit data or string with
%                       path to .mat file with saved fit data.
%
% OUTPUT:
%                dU --> (n x 3) array where n number of T_ref evaluation
%                       points. 
%                       1st column: T_ref
%                       2nd column: dU from measurement
%                       3rd column: dU from optimization

switch nargin
    case 1  % compute dU from results of a finished fit

        if isstruct(varargin{1})
            fit_data = varargin{1};
        elseif ischar(varargin{1})
            fit_data = load(varargin{1});
        end

        % abbreviations
        dsc_data = fit_data.measurements.dsc_data_struct.data;
        index_T_dsc = fit_data.measurements.index_T_dsc;

        p_sim = fit_data.simulation;
        heat_rate = p_sim.heat_rate / 60; % [K/min] -> [K/s]
        T_0 = p_sim.T_0;
        T_end = p_sim.T_end;
        lambda_const = p_sim.lambda_test_setup(1);
        rho_const = p_sim.rho_test_setup(1);
        c_p_const = p_sim.c_p_test_setup(1);
        L1 = p_sim.L1;
            
        p_optim = fit_data.optimization.param_end;

        % compute analytical solution for T_ref
        dt = 0.05 / heat_rate; % fct evaluation every 0.05K
        t = 0:dt:1/heat_rate*(T_end - T_0);
        n = 100;
        a = lambda_const / (c_p_const * rho_const);
        T_ref = analytical_sol(L1,t,n,T_0, heat_rate, a); 

        % get measurement values
        U_dsc = [dsc_data(index_T_dsc(1):index_T_dsc(2),1), dsc_data(index_T_dsc(1):index_T_dsc(2),3)];

        if fit_data.measurements.revMassNorm
            U_dsc(:,2) = U_dsc(:,2) * fit_data.measurements.dsc_data_struct.mass;
        end

    case 4  % compute dU as a sub-step in optimization process to get the 
            % residuum of measured and simulated potential difference.

        p_sim = varargin{1};
        T_ref = varargin{2};
        U_dsc = varargin{3};
        p_optim = varargin{4};
    
    otherwise
    error(['compute_dU needs either 1 argument (finished fit results) or ' ...
           '4 arguments when used in optimization process to get the residuum!'])
end


T_pcm = simulate_1d(p_sim.eval_c_p, p_sim.eval_dc_p, p_sim);

dT = T_ref(:) - T_pcm(:,p_sim.N1);

k = p_sim.get_param_k(p_optim);
dU = polyval(k, dT);

% Note: In the first few seconds T_ref(:) stays constant till the heat of 
% the oven reaches it. Interpolation can't handle non-unique domain of
% definition.
index_T_p5 = find(T_ref(:) > p_sim.T_0 + 5, 1, 'first');

dU_interp = interp1(T_ref(index_T_p5:end), dU(index_T_p5:end), ...
                    U_dsc(:,1), 'linear');
                
% Note: Interpolation necessary because T_ref doesnt increase linearly, at
% least at the beginning, i.e. we dont start at T_0, there's a delay of the
% Constantan.

                
switch nargin
    case 1
        dU_output = [U_dsc(:,1), U_dsc(:,2), dU_interp];

    case 4
        dU_output = dU_interp;

end

