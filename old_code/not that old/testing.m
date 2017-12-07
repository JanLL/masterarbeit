% forward integration with NURBS

L1 = 15.;
L2 = 0.;
L3 = 0.5;

N3 = 50;

T_0 = 10;
T_end = 200;

heat_rate = 10.; % K/min

lambda_test_setup = [23*1, 35.6000, 0.9600];


% c_p parametrization of sample
nrb_order = 4; % nrb_order = 4 equates to C^2

cntrl_pts = [10, 30, 60, 90, 120, 125, 130, 131, 135, 150, 160, 200; ...
             1, 1,  1.1, 1.15, 1.2, 10, 10, 10, 10, 1.52, 1.53, 1.58];
         
% "best" 10K/min fit
cntrl_pts = [0,30,60,90,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,160,180,200; ...
             1.79858157644537,1.73123570378148,1.88939296848618,2.38582410198609,2.83324784260016,2.42651409950467,4.16171664906883,1.13860062268564,5.04919663138140,2.07976511011705,4.83068638838580,3.63392796870435,4.55952727527380,7.61255409096836,2.54742792701244,38.5762115714479,41.9573364402785,11.0500960538204,7.73154658124738,2.65644380692096,2.24473765106880,1.50806399242812,1.38297861923093,1.19283792707801,1.18634381049379,1.21363453226050,1.33993164075706,1.62786023313508,1.76526003655576,2.19534568189373,2.56165109798611,9.17894802559609,2];
         
num_cntrl_pts = size(cntrl_pts,2);

% equidistant knots
knots = [zeros(1,nrb_order), ...
         (1:num_cntrl_pts-nrb_order)/(num_cntrl_pts-nrb_order+1), ...
         ones(1,nrb_order)];

nurbs = nrbmak(cntrl_pts, knots);
dnurbs = nrbderiv(nurbs);                 

tt = 0:0.001:1;
[C, dC] = nrbdeval(nurbs, dnurbs, tt);
     
figure(1)
clf;
%plot(C(1,:), C(2,:))
plot(tt, C(1,:)); hold on   
%plot(tt, dC(1,:))


%%% Newton Method
t0 = 0.5;
Tx = 110.;

h = 0.1;

tk = 0.5;
err = inf;

i = 0;  % iteration counter
tic;
while err > 0.5
    
    [Tk, dTk] = nrbdeval(nurbs, dnurbs, tk); 
    dt = -h*(Tk(1,1)-Tx)/dTk(1,1);
    tk = tk + dt;


    Tknew = nrbeval(nurbs, tk);
    Tknew = Tknew(1,1);
    err = abs(Tknew - Tx);
    %fprintf('%g\t%e\n', tk, err);
    
    i = i+1;

end
toc
fprintf('%d\n', i);
    

return


c_p_sample = {'NURBS', [size(cntrl_pts,2), length(knots)]};

common_args = {'L1', L1, 'L2', L2, 'L3', L3, 'N3', N3, 'T_0', T_0, ...
               'T_end', T_end, 'heat_rate', heat_rate, ...
               'lambda_test_setup', lambda_test_setup, ...
               'c_p_sample', c_p_sample};
p_sim = get_param_sim(common_args{:});

p_optim = [cntrl_pts(1,:), cntrl_pts(2,:), knots];

p_sim = update_c_p(p_sim, p_optim);



% load parameters from fit result
fit_data = load('/home/argo/masterarbeit/fits_data/2017-08-13_15:49:01_407_10Kmin_lsqnonlin/fit_data.mat');
p_sim = fit_data.simulation;


% tic;
T_pcm_1 = simulate_1d(p_sim);
% toc



% Analytical reference solution
lambda_const = p_sim.lambda_test_setup(1);
rho_const = p_sim.rho_test_setup(1);
c_p_const = p_sim.c_p_test_setup(1);

heat_rate_s = p_sim.heat_rate / 60;
dt = 0.05 / heat_rate_s; % fct evaluation every 0.05K
t = 0:dt:1/heat_rate_s*(p_sim.T_end - p_sim.T_0);
n = 100;
a = lambda_const / (c_p_const * rho_const);
x = [0, p_sim.L1];
T_ref = analytical_sol(x,t,n,T_0, heat_rate_s, a); 


fprintf('T_pcm(T_ref = 160 degC) = %g degC\n', T_pcm(find(T_ref > 160, 1, 'first'), p_sim.N1));



