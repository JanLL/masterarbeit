T0 = 30.;
beta = 10. / 60;
L = 25.;

c_p = 0.41;
rho = 8.9;
lambda = 23.;

a = lambda / (c_p * rho);

lambda_sqrt = @(n) ((2*n-1)*pi/(2*L));

u = @(x,t,n) T0 + beta*t(:,:,1) - sum(4.*beta./(pi.*(2.*n-1)) .* 1./(a.*lambda_sqrt(n).^2) .* ...
          (1 - exp(-a.*lambda_sqrt(n).^2.*t)) .* sin(lambda_sqrt(n).*x),3);
       
x = 0:0.05:25;
t = 0:1:300;



N = 1:20;
N = reshape(N,[1,1,length(N)]);
N = repmat(N,[length(t), length(x), 1]);

[X, T] = meshgrid(x,t);
X = repmat(X, [1,1,size(N,3)]);
T = repmat(T, [1,1,size(N,3)]);

sol20 = u(X,T,N);


%  image(u(X,T)','CDataMapping','scaled');
%  colorbar;
%  xlabel('time');
%  ylabel('x');





% Numerical Solution

% simulation data
L1 = 25.;
L2 = 0.;
L3 = 1.;
N3 = 50;

T_0 = 30;
T_end = 200;

heat_rate = 10.; % K/min

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate};
p_sim = get_param_sim(common_args{:});

c_p_params = [144.0009 - 15., ...
                    4.1036 * 5., ...
                    0.0039 + 0.1, ...
                    1.4217 * 0., ...
                    0.0078, ...
                    1.5325];
                
eval_c_p_expl = @(T) p_sim(1).eval_c_p(T, c_p_params);
eval_dc_p_expl = @(T) p_sim(1).eval_dc_p(T, c_p_params);

T_ref = simulate_1d(eval_c_p_expl, eval_dc_p_expl, p_sim(2));


plot(sol20(:,1), sol20(:,1) - sol20(:,end), '--'); hold on
plot(T_ref(:,1), T_ref(:,1) - T_ref(:,end), '--');

xlabel('T_oven [degC]');

