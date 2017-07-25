function [q_pcm_in_interp] = compute_q_pcm_in(varargin)
% TODO: function description, ...


switch nargin
    case 4 % compute q_pcm_in as a sub-step in optimization process to get 
            % the residuum of measured and simulated heat flux into the pcm.
        
    p_sim = varargin{1};
    T_ref = varargin{2};
    q_dsc = varargin{3};
    m_pcm = varargin{4};

    otherwise
    error('Currently just optimization subroutine case with 3 input variables implemented!')

end

% abbreviations
N1 = p_sim.N1;
N3 = p_sim.N3;
L3 = p_sim.L3;
lambda_pcm = p_sim.lambda_test_setup(3);
dt = 0.05 / p_sim.heat_rate * 60; % function evaluation every 0.05K at oven

T_pcm = simulate_1d(p_sim);


q_N1p1_out = - (lambda_pcm * m_pcm)./(N3 * rho_formula(T_pcm(1:end-1,N1+1))) ...
              .*(T_pcm(1:end-1,N1+2) - T_pcm(1:end-1,N1+1)) / (L3/N3)^2;

dQdt = p_sim.eval_c_p(T_pcm(1:end-1,N1+1))*m_pcm/N3 ...
       .* (T_pcm(1:end-1,N1+1) - T_pcm(2:end,N1+1)) / dt;
   

q_pcm_in = dQdt + q_N1p1_out;

% Note: In the first few seconds T_ref(:) stays constant till the heat of 
% the oven reaches it. Interpolation can't handle non-unique domain of
% definition.
index_T_p5 = find(T_ref(:) > p_sim.T_0 + 5, 1, 'first');
q_pcm_in_interp = interp1(T_ref(index_T_p5:end-1), q_pcm_in(index_T_p5:end), ...
                    q_dsc(:,1), 'linear');

end

