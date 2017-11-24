% %% Compare c_p fits for heat rate 0,3 K/min for different L1/L3/N1/N3
% 
% fit_path_root_list = {'2017-10-18_12:11:40_407_L1=40_L3=0.01_N1=200_N3=50', ...
%                       '2017-10-16_04:45:00_407_L1=40_L3=0.1_N1=200_N3=50', ...
%                       '2017-10-18_14:41:19_407_L1=40_L3=0.2_N1=200_N3=50', ...
%                       '2017-10-18_13:28:22_407_L1=5_L3=0.01_N1=200_N3=50', ...
%                       '2017-10-18_10:59:34_407_L1=80_L3=0.1_N1=200_N3=50'};
% 
%                   
% T_domain = 100:0.005:160;
% figure(); hold on
% 
% for j=1:length(fit_path_root_list)
% 
%     fit_path_root = strcat('/home/argo/masterarbeit/fits_data/', fit_path_root_list{j}, '/');
% 
%     
%     file_list = dir(fit_path_root);
% 
%     isub = [file_list(:).isdir]; %# returns logical vector
%     nameSubDirs = {file_list(isub).name}';
%     nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];
% 
% 
%     for i=1:length(nameSubDirs)
% 
%         if contains(nameSubDirs{i}, '0,3Kmin')
% 
%             fit_data_path = strcat(fit_path_root, nameSubDirs{i}, '/fit_data.mat');
%             
%             fit_data = open(fit_data_path);
% 
%             c_p = c_p_gauss_linear_comb(T_domain, fit_data.optimization.p_optim_end);
%             plot(T_domain, c_p, 'DisplayName', strcat(num2str(fit_data.measurement.dsc_data.Tinfo.Tstep), ' K/min'))
% 
%         end
% 
%     end
% 
% end
% 
% 
% return


%% Plots c_p and heat flux from fit in an appropriate way

c_p_fig = open('c_p(T).fig');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
children = get(gca, 'Children');
children.LineWidth = 2;
%set(gca,'FontSize',24)
set(gca,'xlim', [20 max(children.XData)]);
set(gca,'xlim', [110 145]);
ylabel('c_p [mJ/(mg*K)]');
title('');

print(c_p_fig, 'c_p(T)', '-dpng', '-r200');
close;
%%
c_p_fig = open('q_pcm_in(T_ref).fig');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
children = get(gca, 'Children');
children(1).LineWidth = 2.;
children(2).LineWidth = 2.;
children(3).LineWidth = 2.;
children(2).LineStyle = '--';
set(gca,'FontSize',24)
set(gca,'xlim', [min(children(1).XData) - mod(min(children(1).XData),20) 160]);  %
title('');
ylabel('\Phi_q^{PCM,in}[mW]')


print(c_p_fig, 'q_pcm_in(T_ref)', '-dpng', '-r200');
close;


return


%% Plot heat flux measurements and c_p computed from DIN formula

dsc_list = DSC204_readFiles(['/home/argo/masterarbeit/', ...
    'DSC204_F1_Phoenix_Messungen/Messungen/Messungen/', ...
    'ExpDat_16-407-3_mitKorr_*Kmin_H.csv']);

fig1 = figure(1); hold on
ax1 = gca();
fig2 = figure(2); hold on
ax2 = gca();

for i=1:length(dsc_list)

    dsc = dsc_list(i);
    q_meas = dsc.data(:,3) ./ dsc.data(:,4) * dsc.mass;
    legend_str = [num2str(dsc.Tinfo.Tstep), ' K/min'];
    plot(ax1, dsc.data(:,1), q_meas, 'DisplayName', legend_str, ...
        'LineWidth', 2);
    
    c_p = calc_cp(dsc);
    legend_str = [num2str(dsc.Tinfo.Tstep), ' K/min'];
    plot(ax2, c_p(:,1), c_p(:,2), 'DisplayName', legend_str, ...
        'LineWidth', 2);
    
    
end

figure(1);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
set(gca,'FontSize',24);
xlabel('T_{ref}[degC]');
ylabel('\eta^{\Phi}[mW]');
set(gca, 'xlim', [30 160]);
legend('show', 'location', 'northwest');

print(fig1, 'heat_flux_measurement', '-dpng', '-r200');
close();


figure(2);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
set(gca,'FontSize',24);
xlabel('T [degC]');
ylabel('c_p [mJ/(mg*K)]');
set(gca, 'xlim', [30 160]);
legend('show', 'location', 'northwest');

print(fig2, 'c_p_DIN_formula', '-dpng', '-r200');
close();


return




%% Plot analytical solution of reference temperature
%  (wird in Vortrag/Thesis genutzt als Illustration wie man von
%  Referenztemperaturen mittels Newton-Verfahren auf Messzeitpunkte kommt.)
L1 = 20.;
t = 0:0.1:1500;
n = 1;
T_0 = 10.;
heat_rate_s = 10/60;
a_Const = 1.;

T_ref = analytical_sol(L1,t,n,T_0, heat_rate_s, a_Const);
T_furnace = T_0 + heat_rate_s * t;

figure(1); hold on

plot(t, T_furnace, 'DisplayName', 'T_{furnace}', 'linewidth', 2.);
plot(t, T_ref, 'DisplayName', 'T_{ref}', 'linewidth', 2.)
xlabel('time')
ylabel('Temp')
legend('show', 'location', 'northwest')





%% Convergence Rate Plots

fit_data = load('fit_data.mat');


%  || x_i - x^* ||^2
n_optim_process = size(fit_data.optimization.p_optim_process,1);
p_optim_process_norm1 = zeros(n_optim_process,1);
for i=1:n_optim_process
    
    p_optim_process_norm1(i) = sum((fit_data.optimization.p_optim_process(i,:) - fit_data.optimization.p_optim_end).^2);
    
end

figure();
plot(p_optim_process_norm1, 'x'); hold on

c = zeros(n_optim_process-1,1);
for i=1:n_optim_process-1
   
    c(i) = p_optim_process_norm1(i+1) / p_optim_process_norm1(i);
    
end

figure();
plot(c, 'x')






% Testing of fraser-suzuki peak
h  =  30.0;
r  =  30.0;
wr =  15.0;
sr =   0.3;
z  = 125.0;
b  =   2.0;

border = z - wr*sr/(sr^2 - 1)

c_p = @(T) 30*exp(-log(30)/log(0.3)^2 * (log(1+(T-125)*(0.3^2-1)/(15*0.3))).^2);

syms T;
dc_p = diff(c_p(T), T);

T = (border-0.001):0.0001:(border-0.0001);
plot(T, c_p(T))

return

% testing of old own c_p parametrization with arctan

syms T;
p = sym('p', [6, 1]);
dc_p_formula = matlabFunction(diff(c_p_formula(T, p), T), 'Vars', [T;p]);
dc_p_formula = @(T,p) dc_p_formula(T,p(1),p(2),p(3),p(4),p(5),p(6));

dc_p_derivation = @(T,p) -(p(2)*p(4)*exp(-p(3)*(T-p(1)).^2))./(1+(p(4)*(T-p(1))).^2) ...
                         - 2*p(2)*p(3)*(T-p(1)).*exp(-p(3)*(T-p(1)).^2) ...
                         .* (atan(-p(4)*(T-p(1)))+pi/2) + p(5);


p = [130., 10, 0.01, 0.5, 0.003, 2];

T = 30:0.01:160;
c_p = c_p_formula(T, p);
dc_p = dc_p_formula(T, p);

dc_p_der = dc_p_derivation(T, p);

figure(1)
clf;
plot(T, c_p); hold on
plot(T, dc_p)
plot(T, dc_p_der, '--')

return



%% testing of frazer-suzuki parametrization

r  =   20;
h  =  40;
z  = 120.0;
wr =  25;
sr =   0.4;
x  = 3:0.2:250;


[y, dy] = frazersuzuki(r, h, z, wr, sr, x);


clf;
plot(x, y); hold on
plot(x, dy)



%% Plotting relative error of discretization grid size in a pretty way
fig = open('relError_comparison_hetero_grid_200_2500_Ag.fig');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
children = get(gca, 'Children');
children.LineWidth = 2;
set(gca,'FontSize',24);
xlabel('measurement points');
ylabel('Abs. Relative Error of \Phi_q^{PCM, in}')
print(fig, 'relErr_discretization_grid_200_2500', '-dpng', '-r200');
close();


%% Plotting computed c_p and corresponding integral, the enthalpy

enthalpy_data = load('enthalpy_test.txt');
T_tilde = 30.;

T = enthalpy_data(:,1);
c_p = enthalpy_data(:,2);
h = enthalpy_data(:,3);

params = [16.6062086431708,	1.92927098591158,	128.850918977868, ...
		34.1327529190694,	1.27007399531043,	130.601505840094, ...
		  4.66622071658818,	3.08315397789580,	126.387021606664, ...
		  0.330799303747586,	6.28205711993873,	138.931914796423, ...	
		 -1.99975880363326,	1.00794170224404,	127.886002738954, ...
		2.00575452880237,	5.23818264316667,	122.475803758605, ...
		8.69795539050226,	0.764167761467260,	131.898468412807, ...
			0.,				26.4978102269972,	286.560297948564, ...
			0.,				11.4117833489350,	111.924264531493, ...
            0.,				35.7646524748368,	95.2324216508870, ...
			0.,				0.];

c_p1 = c_p_gauss_linear_comb(T, params);
        
h1 = 0.;
for i=1:3:30
    
    A_i = params(i);
    sigma_i = params(i+1);
    T_i = params(i+2);
    
    h1 = h1 + sqrt(pi)/2. * A_i * sigma_i ...
        * (erf((T - T_i)/sigma_i) - erf((T_tilde - T_i)/sigma_i));
    
end

h1 = h1 + 0.01*params(31) / 2. * (T.^2 - T_tilde^2) + params(32) * (T - T_tilde);


figure(5); hold on
plot(T, c_p);
plot(T, c_p1, '--');

figure(6); hold on
plot(T, h);
plot(T, h1, '--');


%% Variable grid with logistic growth function
N1 = 300;
N3 = 50;
L1 = 40.;
L3 = 0.1;

x_bottom = 18/20*N1;
x_top = 19/20*N1;
threshold = 0.95;
% threshold = 0.5; % equidistant grid for all N1, N3

N = N1+N3;
x0 = 1/2*(x_top + x_bottom);
gamma = log(1/threshold - 1) / (1/2*(x_bottom - x_top));
dx_pcm = L3 / (N3-1);
%dx_pcm = L1 / (N-1);


dx_Ag = (L1 + L3 - (N-1)*dx_pcm) ./ sum(1./(1+exp(gamma*((0:N-2) - x0)))) + dx_pcm;


f = @(x) (dx_Ag - dx_pcm)./(1 + exp(gamma*(x-x0))) + dx_pcm;
%F = @(x) (dx_Ag - dx_pcm) * (x + 1./gamma*log((exp(-gamma*x0)+1)./(exp(gamma*(x-x0))+1))) ...
%         + x*dx_pcm;

x = 0:1:N-1;

figure(1); clf;
plot(x,f(x), 'x')
xlabel('discretization point')
ylabel('Grid size \Delta x')
set(gca,'FontSize',12);


figure(2); clf;
plot(x,cumsum(f(x)), 'x')
xlabel('discretization point')
ylabel('Physical Grid $$\tilde{x}$$', 'Interpreter','latex')
set(gca,'FontSize',12);

sum(f(0:N-2))
sum(f(0:N1-2))



%% New variable grid with equidistant Rechengitter [0,1]
% falsche Laengen bei N1 und N1+N3 am Ende weil Integral benutzt statt
% Summe.

% Parameters to set
L1 = 40.;
L3 = 1.;

N = 300;
n_pcm = 0.2;
n_tr = 0.7;
n_m = 0.01;
t = 0.99;

b = 1 - n_pcm - n_m - n_tr/2;
gamma = 2/n_tr * log(t/(1-t));

Q11 = 1 + 1/gamma * log((exp(-gamma*b) + 1) / (exp(gamma*(1-b)) + 1));
Q12 = -1/gamma * log((exp(-gamma*b) + 1) / (exp(gamma*(1-b)) + 1));
Q21 = 1 - n_pcm + 1/gamma * log((exp(-gamma*b) + 1) / (exp(gamma*(1-n_pcm-b)) + 1));
Q22 = -1/gamma * log((exp(-gamma*b) + 1) / (exp(gamma*(1-n_pcm-b)) + 1));

Q = [Q11, Q12;
     Q21, Q22];

 dx = Q\[L1+L3;L1];
 dx_Ag = dx(1);
 dx_pcm = dx(2);
 
chi = @(x_tilde) (dx_Ag - dx_pcm)*(x_tilde + 1/gamma * log((exp(-gamma*b)+1) ./ (exp(gamma*(x_tilde-b))+1))) + dx_pcm*x_tilde; 
dchi = @(x_tilde) (dx_Ag - dx_pcm) ./ (exp(gamma*(x_tilde-b)) + 1) + dx_pcm;
 
x_tilde_domain = linspace(0,1,N);

figure(1)
plot(x_tilde_domain*N, chi(x_tilde_domain), 'x')

figure(2); clf; hold on;
%plot(x_tilde_domain*N, gradient(chi(x_tilde_domain), x_tilde_domain*N), 'x')
plot(N*x_tilde_domain, dchi(x_tilde_domain), 'x')



%% New variable grid with equidistant Rechengitter {0,1,...,N-1}

% Parameters to set
L1 = 40.;
L3 = 0.1;

N1 = 300;
N3 = 50;
n_pcm = N3 / (N1 + N3);
n_tr = 0.3;
n_m = 0.01;
t = 0.999;

N = N1 + N3;
N_pcm = N * n_pcm;
N_tr = N * n_tr;
N_m = N * n_m;

b = N-1 - N_pcm - N_m - N_tr/2;
gamma = 2/N_tr * log(t/(1-t));

W11 = sum(1./(exp(gamma*((0:N1-2) - b)) + 1));
W12 = N1 - 1 - sum(1./(exp(gamma*((0:N1-2) - b)) + 1));
W21 = sum(1./(exp(gamma*((0:N1+N3-2) - b)) + 1));
W22 = N1 + N3 - 1 - sum(1./(exp(gamma*((0:N1+N3-2) - b)) + 1));

W = [W11, W12;
     W21, W22];

 dx = W\[L1;L1+L3];
 dx_Ag = dx(1);
 dx_pcm = dx(2);
 
chi = @(x_tilde) (dx_Ag - dx_pcm)*(x_tilde + 1/gamma * log((exp(-gamma*b)+1) ./ (exp(gamma*(x_tilde-b))+1))) + dx_pcm*x_tilde; 
dchi = @(x_tilde) (dx_Ag - dx_pcm) ./ (exp(gamma*(x_tilde-b)) + 1) + dx_pcm;
 
x_tilde_domain = 0:N-1;

%figure(3); clf; hold on
%plot(x_tilde_domain, cumsum(dchi(x_tilde_domain)), 'x')

figure(4); clf; hold on;
plot(x_tilde_domain, dchi(x_tilde_domain), 'x')


sum(dchi(0:N1-2))
sum(dchi(0:N1+N3-2))



%% Fraser-Suzuki Beispielbild

% c_p parametrization with Fraser-Suzuki-Peak
h  =  10.0;
r  =  2.0;
wr =  5.0;
sr =   0.3;
z  = 130.0;
m  = 0.003;
b  =   2.0;
p_fraser_suzuki = [h, r, wr, sr, z, m, b].';


T = 30:0.001:160;
c_p = c_p_fs(T, p_fraser_suzuki);
fig = figure(5); cla;
set(gca,'FontSize',24);
plot(T,c_p, 'linewidth', 1);
xlabel('T [degC]');
ylabel('c_p [mJ/(mg*K)]')


print(fig, '/home/argo/masterarbeit/thesis/images/fraser_suzuki_example', '-dpng', '-r200');
close();


%% Plot old NURBS results in a pretty way

% Jacobian
% fit_data = load('fit_data.mat');
% 
% fig = open('jac_output.fig');
% set(gca,'FontSize',12);
% 
% print(fig, '/home/argo/masterarbeit/thesis/images/NURBS_jac', '-dpng', '-r200');
% 
% % Zoom in interesting part
% xlim([0, 30])
% ylim([124, 132])
% 
% print(fig, '/home/argo/masterarbeit/thesis/images/NURBS_jac_zoom', '-dpng', '-r200');
% close();


% c_p(T)
% fig = open('c_p(T).fig');
% 
% set(gca,'FontSize',12);
% 
% children = get(gca,'Children');
% delete(children(1)) % delete c_p calc from formula
% 
% ylabel('c_p [mJ/(mg*K)]')
% 
% print(fig, '/home/argo/masterarbeit/thesis/images/NURBS_c_p(T)', '-dpng', '-r200');


% heat flux
fig = open('q_pcm_in(T_ref).fig');

children = get(gca, 'Children');
children(1).LineWidth = 1.;
children(2).LineWidth = 1.;
children(3).LineWidth = 1.;
children(2).LineStyle = '--';

set(gca,'FontSize',12);

ylabel('\Phi_q^{PCM,in} [mW]')

legend('off')
legend('show', 'location', 'northwest');

print(fig, '/home/argo/masterarbeit/thesis/images/NURBS_heat_flux', '-dpng', '-r200');
close;


%% Plot smearing effect in simulation in a pretty way

% heat fluxes
fig = open('heat_flux_smearing_simulation.fig')
set(gca,'FontSize',12);

children = get(gca, 'Children');
children(1).LineWidth = 1.;
children(2).LineWidth = 1.;
children(3).LineWidth = 1.;
children(4).LineWidth = 1.;
children(5).LineWidth = 1.;
children(6).LineWidth = 1.;
children(7).LineWidth = 1.;

print(fig, '/home/argo/masterarbeit/thesis/images/smearing_effect_simulation_heat_fluxes', '-dpng', '-r200');
close();

% c_p(T)
fig = open('c_p(T).fig');

set(gca,'FontSize',12);

children = get(gca, 'Children');
children(1).LineWidth = 1.;

ylabel('c_p [mJ/(mg*K)]')

print(fig, '/home/argo/masterarbeit/thesis/images/smearing_effect_simulation_c_p', '-dpng', '-r200');

close();



%% Plot grid error with equidistant grid as reference
%  modify n_tr

root_dir = '/home/argo/masterarbeit/simulationen-data/grid_error/';

T_ref_file = [root_dir, 'Temp_equidistant_grid_N1=20000_N3=50_L1=40,00_L3=0,10.mat'];
T_ref_data = load(T_ref_file);
T_ref = T_ref_data.T(:,20000);

T_test_dirs = {'Temp_n_pcm=0,14286_n_tr=0_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,3_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,5_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,7_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat'};

figure(20); clf; ax1=gca; hold on;
figure(21); clf; ax2=gca; hold on;

           
for i=1:length(T_test_dirs)
  
    T_test_data = load([root_dir, T_test_dirs{i}]);
    T_test = T_test_data.T(:,T_test_data.N1);
    
    relErr = abs(1 - T_test ./ T_ref);
    
    plot(ax1, relErr, 'DisplayName', sprintf('%s', num2str(T_test_data.n_tr)));
    
    %T_test_data.N1
    
    spatial_grid = build_grid(T_test_data.N1, T_test_data.N3, T_test_data.L1, T_test_data.L3, T_test_data.n_tr, T_test_data.n_m, T_test_data.t);
    
    plot(ax2, spatial_grid, 'x', 'DisplayName', sprintf('%s', num2str(T_test_data.n_tr)));
    
end

legend(ax1, 'show', 'location', 'northwest')
legend(ax2, 'show', 'location', 'northeast')


%%
%  modify N1

root_dir = '/home/argo/masterarbeit/simulationen-data/grid_error/';

T_ref_file = [root_dir, 'Temp_equidistant_grid_N1=20000_N3=50_L1=40,00_L3=0,10.mat'];
T_ref_data = load(T_ref_file);
T_ref = T_ref_data.T(:,20000);

T_test_dirs = {'Temp_n_pcm=0,33333_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=100_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,090909_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=500_N3=50.mat', ...
               'Temp_n_pcm=0,047619_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=1000_N3=50.mat'};

figure(20); clf; ax1=gca; hold on;
figure(21); clf; ax2=gca; hold on;

           
for i=1:length(T_test_dirs)
  
    T_test_data = load([root_dir, T_test_dirs{i}]);
    T_test = T_test_data.T(:,T_test_data.N1);
    
    relErr = abs(1 - T_test ./ T_ref);
    
    plot(ax1, relErr, 'DisplayName', sprintf('%s', num2str(T_test_data.N1)));
    
    %T_test_data.N1
    
    spatial_grid = build_grid(T_test_data.N1, T_test_data.N3, T_test_data.L1, T_test_data.L3, T_test_data.n_tr, T_test_data.n_m, T_test_data.t);
    
    plot(ax2, spatial_grid, 'x', 'DisplayName', sprintf('%s', num2str(T_test_data.N1)));
    
end

legend(ax1, 'show', 'location', 'northwest')
legend(ax2, 'show', 'location', 'northeast')




%%
%  modify N

root_dir = '/home/argo/masterarbeit/simulationen-data/grid_error/';

T_ref_file = [root_dir, 'Temp_equidistant_grid_N1=20000_N3=50_L1=40,00_L3=0,10.mat'];
T_ref_data = load(T_ref_file);
T_ref = T_ref_data.T(:,20000);

T_test_dirs = {'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=160_N3=40.mat', ...
               'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=240_N3=60.mat', ...
               'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=400_N3=100.mat', ...
               'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=800_N3=200.mat'};

figure(20); clf; ax1=gca; hold on;
figure(21); clf; ax2=gca; hold on;

           
for i=1:length(T_test_dirs)
  
    T_test_data = load([root_dir, T_test_dirs{i}]);
    T_test = T_test_data.T(:,T_test_data.N1);
    
    relErr = abs(1 - T_test ./ T_ref);
    
    plot(ax1, relErr, 'DisplayName', sprintf('%s', num2str(T_test_data.N1 + T_test_data.N3)));
    
    spatial_grid = build_grid(T_test_data.N1, T_test_data.N3, T_test_data.L1, T_test_data.L3, T_test_data.n_tr, T_test_data.n_m, T_test_data.t);
    
    plot(ax2, spatial_grid, 'x', 'DisplayName', sprintf('%s', num2str(T_test_data.N1 + T_test_data.N3)));
    
end

legend(ax1, 'show', 'location', 'northwest')
legend(ax2, 'show', 'location', 'northeast')





