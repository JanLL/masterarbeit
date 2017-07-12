% normal forward integration

L1 = 25.;
L2 = 0.;
L3 = 1.;

N1 = 5000;
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


% Numerical sample solution
tic;
T_pcm = simulate_1d(p_sim.eval_c_p, p_sim.eval_dc_p, p_sim);
toc


% Analytical reference solution
lambda_const = p_sim.lambda_test_setup(1);
rho_const = p_sim.rho_test_setup(1);
c_p_const = p_sim.c_p_test_setup(1);

heat_rate_s = heat_rate / 60;
dt = 0.05 / heat_rate_s; % fct evaluation every 0.05K
t = 0:dt:1/heat_rate_s*(T_end - T_0);
n = 100;
a = lambda_const / (c_p_const * rho_const);
T_ref = analytical_sol(L1,t,n,T_0, heat_rate_s, a); 

N1 = p_sim.N1;
N2 = p_sim.N2;
N3 = p_sim.N3;

dx = ones(N1+N2+N3, 1);
dx(1:N1)        = L1/N1;
dx(N1+1:N1+N2)  = L2/N2;
dx(N1+N2+1:end) = L3/N3;
pcm_const = p_sim.lambda_test_setup(3);


plot(T_ref,T_ref-T_pcm(:,N1)); hold on



% dx unkritisch solange ueberall gleich, ansonsten wirds etwas
% komplizierter...
%q_const = -lambda_const * (T_pcm(:,N1) - T_pcm(:,N1-1))/dx(N1-1);
%q_pcm = pcm_const * (T_pcm(:,N1+2) - T_pcm(:,N1+1))/dx(N1+1);
%q_pcm = -pcm_const * sum(diff(T_pcm(:,N1+1:end),1,2) / dx(N1+1),2);

%plot(T_ref, q_const, 'DisplayName', 'Constantan'); hold on
%plot(T_ref, q_pcm, 'DisplayName', 'PCM')
%legend('show', 'location', 'northwest')

% Kleine Abweichung, denke aufgrund Waermestromanteil beim Uebergang
% Constantan -> PCM, wo allerdings der Waermeleitkoeff. nicht bekannt ist.
% Man koennte mit der Abweichung und dem geg. Temp.-grad. diesen
% ausrechnen...


