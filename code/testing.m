% convergence of spatial discretization N1

L1 = 25.;
L2 = 0.;
L3 = 1.;

N1 = 6;
N2 = 0;
N3 = 50;

T_0 = 10;
T_end = 300;

heat_rate = 10.; % K/min

lambda_test_setup = [23*1, 35.6000, 0.9600];


common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N1', N1, 'N2', N2, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup};
p_sim = get_param_sim(common_args{:});


c_p_params = [144.0009, ...
                    4.1036, ...
                    0.0039, ...
                    1.4217, ...
                    0.0078, ...
                    1.5325];

p_sim = update_c_p(p_sim, c_p_params);


T_pcm_old = simulate_1d(p_sim.eval_c_p, p_sim.eval_dc_p, p_sim);

error = 1;
N1_list = [N1];

i = 1;
for i=1:50

    N1_old = N1;
    N1 = N1 + 2;
    p_sim.N1 = N1;
    
    T_pcm_new = simulate_1d(p_sim.eval_c_p, p_sim.eval_dc_p, p_sim);
    
    error = max(abs((T_pcm_new(:,N1) - T_pcm_old(:,N1_old)) ./ T_pcm_new(:,N1)));
    
    error_list(i) = error;
    %i = i+1;
    
    T_pcm_old = T_pcm_new;
    
end


semilogy((1:length(error_list))*2 + 6, error_list, 'x')
xlabel('N1');
ylabel('max\{|1 - T_{pcm}(N_1) / T_{pcm}(N_1^+)|\}');
title(sprintf('Convergence rate of increasing Constantan\ndiscretization number'));













