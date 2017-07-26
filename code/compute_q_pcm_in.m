function [q_pcm_in_output] = compute_q_pcm_in(varargin)
% TODO: function description, ...


switch nargin
    case 1 % compute q_pcm_in from results of a finished fit
        
        if isstruct(varargin{1})
            fit_data = varargin{1};
        elseif ischar(varargin{1})
            fit_data = load(varargin{1});
        end

        % abbreviations
        dsc = fit_data.measurements.dsc_data_struct;
        index_T_dsc = fit_data.measurements.index_T_dsc;
        m_pcm = dsc.mass;

        p_sim = fit_data.simulation;
        heat_rate = p_sim.heat_rate / 60; % [K/min] -> [K/s]
        T_0 = p_sim.T_0;
        T_end = p_sim.T_end;
        lambda_const = p_sim.lambda_test_setup(1);
        rho_const = p_sim.rho_test_setup(1);
        c_p_const = p_sim.c_p_test_setup(1);
        L1 = p_sim.L1;
        
        optim_type_int = str2double(fit_data.optimization.optim_type(end));
            
        % compute analytical solution for T_ref
        dt = 0.05 / heat_rate; % fct evaluation every 0.05K
        t = 0:dt:1/heat_rate*(T_end - T_0);
        n = 100;
        a = lambda_const / (c_p_const * rho_const);
        T_ref = analytical_sol(L1,t,n,T_0, heat_rate, a); 

        % get measurement values
        q_dsc = [dsc.data(index_T_dsc(1):index_T_dsc(2),1), ...
                    dsc.data(index_T_dsc(1):index_T_dsc(2),3) ...
                 ./ dsc.data(index_T_dsc(1):index_T_dsc(2),4)];
        
        if fit_data.measurements.revMassNorm
            q_dsc(:,2) = q_dsc(:,2) * fit_data.measurements.dsc_data_struct.mass;
        end
        
    
    case 5 % compute q_pcm_in as a sub-step in optimization process to get 
            % the residuum of measured and simulated heat flux into the pcm.

        p_sim = varargin{1};
        T_ref = varargin{2};
        q_dsc = varargin{3};
        m_pcm = varargin{4};
        optim_type_int = varargin{5};

    otherwise
    error('Currently just optimization subroutine case with 5 input variables implemented!')

end

% abbreviations
N1 = p_sim.N1;
N3 = p_sim.N3;
L3 = p_sim.L3;
lambda_pcm = p_sim.lambda_test_setup(3);
dt = 0.05 / p_sim.heat_rate * 60; % function evaluation every 0.05K at oven

T_pcm = simulate_1d(p_sim);

if optim_type_int == 1
    % Computation of q_pcm_in via combination of heat change aund fourier's law
    q_N1p1_out = - (lambda_pcm * m_pcm)./(N3 * rho_formula(T_pcm(1:end-1,N1+1))) ...
                  .*(T_pcm(1:end-1,N1+2) - T_pcm(1:end-1,N1+1)) / (L3/N3)^2;

    dQdt = p_sim.eval_c_p(T_pcm(1:end-1,N1+1))*m_pcm/N3 ...
           .* (T_pcm(2:end,N1+1) - T_pcm(1:end-1,N1+1)) / dt;

    q_pcm_in = dQdt + q_N1p1_out;

    % figure(4)
    % plot(q_N1p1_out, 'DisplayName', 'q_{N1p1}^{out}'); hold on
    % plot(dQdt, 'DisplayName', 'dQdt');
    % legend('show');

elseif optim_type_int == 2   
    % Alternative computation of q_pcm_in via sum of all dQdt in PCM
    dQdt = p_sim.eval_c_p(T_pcm(1:end-1,N1+1:end))*m_pcm/N3 ...
           .* (T_pcm(2:end,N1+1:end) - T_pcm(1:end-1,N1+1:end)) / dt;

    q_pcm_in = sum(dQdt,2);
end

% Note: In the first few seconds T_ref(:) stays constant till the heat of 
% the oven reaches it. Interpolation can't handle non-unique domain of
% definition.
index_T_p5 = find(T_ref(:) > p_sim.T_0 + 5, 1, 'first');
q_pcm_in_interp = interp1(T_ref(index_T_p5:end-1), q_pcm_in(index_T_p5:end), ...
                    q_dsc(:,1), 'linear');

                
switch nargin
    case 1
        q_pcm_in_output = [q_dsc(:,1), q_dsc(:,2), q_pcm_in_interp];

    case 5
        q_pcm_in_output = q_pcm_in_interp;
                
                
end

