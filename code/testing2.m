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

%dsc_list = DSC204_readFiles(['/home/argo/masterarbeit/', ...
%    'DSC204_F1_Phoenix_Messungen/Messungen/Messungen/', ...
%    'ExpDat_16-407-3_mitKorr_*Kmin_H.csv']);

dsc_list = {'ExpDat_16-407-3_mitKorr_0,3Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_0,6Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_1,25Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_2,5Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_5Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_10Kmin_H.csv', ...
            'ExpDat_16-407-3_mitKorr_20Kmin_H.csv'};

fig1 = figure(1); hold on
ax1 = gca();
fig2 = figure(2); hold on
ax2 = gca();

for i=1:length(dsc_list)

    dsc = DSC204_readFile(dsc_list{i});
    
    q_meas = dsc.data(:,3) ./ dsc.data(:,4) * dsc.mass;
    legend_str = [num2str(dsc.Tinfo.Tstep), ' K/min'];
    plot(ax1, dsc.data(:,1), q_meas, 'DisplayName', legend_str, ...
        'LineWidth', 1.3);
    
    c_p = calc_cp(dsc);
    legend_str = [num2str(dsc.Tinfo.Tstep), ' K/min'];
    plot(ax2, c_p(:,1), c_p(:,2), 'DisplayName', legend_str, ...
        'LineWidth', 1.3);
    
    
end

save_root_dir = '/home/argo/masterarbeit/thesis/images/'

figure(1);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
set(gca,'FontSize',24);
xlabel('T_{ref}[degC]');
ylabel('\eta^{\Phi}[mW]');
set(gca, 'xlim', [30 160]);
legend('show', 'location', 'northwest');

print(fig1, [save_root_dir, 'heat_flux_measurement'], '-dpng', '-r200');
close();


figure(2);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
set(gca,'FontSize',24);
xlabel('T [degC]');
ylabel('c_p [mJ/(mg*K)]');
set(gca, 'xlim', [30 160]);
legend('show', 'location', 'northwest');

print(fig2, [save_root_dir, 'c_p_DIN_formula'], '-dpng', '-r200');
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

root_dir = '/home/argo/masterarbeit/fits_data/2017-11-19_20:25:11_407_10Kmin_lsqnonlin/'
fit_data = load([root_dir, 'fit_data.mat']);


% Jacobian

fig = open('jac_output.fig');
set(gca,'FontSize',12);

% print(fig, '/home/argo/masterarbeit/thesis/images/NURBS_jac', '-dpng', '-r200');

% Zoom in interesting part
% xlim([0, 30])
% ylim([124, 132])

% print(fig, '/home/argo/masterarbeit/thesis/images/NURBS_jac_zoom', '-dpng', '-r200');
% close();


% c_p(T)
% fig = open('c_p(T).fig');
% 
% set(gca,'FontSize',12);
% 
% children = get(gca,'Children');
% delete(children(1)) % delete c_p calc from formula
% children(2).LineWidth = 1.3;
% 
% ylabel('c_p [mJ/(mg*K)]')
% box off  % entfernt den Rahmen um die figure
% 
% print(fig, '/home/argo/masterarbeit/thesis/images/NURBS_c_p(T)', '-dpng', '-r200');
% close;

% heat flux
% fig = open('q_pcm_in(T_ref).fig');
% 
% children = get(gca, 'Children');
% children(1).LineWidth = 1.3;
% children(2).LineWidth = 1.3;
% children(3).LineWidth = 1.3;
% children(2).LineStyle = '--';
% 
% set(gca,'FontSize',12);
% 
% ylabel('\Phi_q^{PCM,in} [mW]')
% 
% legend('off')
% legend('show', 'location', 'northwest');
% 
% print(fig, '/home/argo/masterarbeit/thesis/images/NURBS_heat_flux', '-dpng', '-r200');
% close;


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

T_ref_file = [root_dir, 'Temp_equidistant_grid_N1=20000_N3=50_L1=40,00_L3=0,10_1e-7.mat'];
% T_ref_file = [root_dir, 'Temp_N1=20000_N3=50_L1=40.00_L3=0.10_relTol=1e-8.mat'];

T_ref_data = load(T_ref_file);
T_ref = T_ref_data.T(:,20000);

T_test_dirs = {'Temp_n_pcm=0,14286_n_tr=0_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,3_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,5_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,7_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat'};           

fig1 = figure(20); clf; ax1=gca; set(ax1, 'YScale', 'log'); hold on;
fig2 = figure(21); clf; ax2=gca; set(ax2, 'YScale', 'log'); hold on;
           
for i=1:length(T_test_dirs)
  
    T_test_data = load([root_dir, T_test_dirs{i}]);
    T_test = T_test_data.T(:,T_test_data.N1);
    
    relErr = abs(1 - T_test ./ T_ref);
    
    plot(ax1, T_test_data.T_ref_dsc, relErr, 'Linewidth', 1.3, 'DisplayName', sprintf('n_{tr}=%s', num2str(T_test_data.n_tr)));
    
    %T_test_data.N1
    
    spatial_grid = build_grid(T_test_data.N1, T_test_data.N3, T_test_data.L1, T_test_data.L3, T_test_data.n_tr, T_test_data.n_m, T_test_data.t);
    
    plot(ax2, spatial_grid, 'x', 'DisplayName', sprintf('n_{tr}=%s', num2str(T_test_data.n_tr)));
    
end

legend(ax1, 'show', 'location', 'northwest', 'Orientation', 'horizontal')
xlabel(ax1, 'T_{ref} [degC]')
ylabel(ax1, '|Relative error|')
set(ax1,'FontSize',12);
ylim(ax1, [2e-7 1e-4]);
set(fig1, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.43, 0.6]); 
print(fig1, '/home/argo/masterarbeit/simulationen-data/grid_error/mod_n_tr_relErr', '-dpng', '-r200');


legend(ax2, 'show', 'location', 'southwest')
xlabel('$$\tilde{x}$$, Computation grid', 'Interpreter', 'latex')
ylabel('$$\frac{\partial \chi}{\partial \tilde{x}}(\tilde{x}) = \Delta x$$ [mm]', 'Interpreter', 'latex')
set(ax2,'FontSize',12);
set(fig2, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.43, 0.6]); 
print(fig2, '/home/argo/masterarbeit/simulationen-data/grid_error/mod_n_tr_gridsize', '-dpng', '-r200');



%%
%  modify N1

root_dir = '/home/argo/masterarbeit/simulationen-data/grid_error/';

T_ref_file = [root_dir, 'Temp_N1=20000_N3=50_L1=40.00_L3=0.10_relTol=1e-7_new.mat'];
% T_ref_file = [root_dir, 'Temp_N1=20000_N3=50_L1=40.00_L3=0.10_relTol=1e-8.mat'];

T_ref_data = load(T_ref_file);
T_ref = T_ref_data.T(:,20050);

T_test_dirs = {'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=200_N3=50.mat', ...
               'Temp_n_pcm=0,14286_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=300_N3=50.mat', ...
               'Temp_n_pcm=0,090909_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=500_N3=50.mat', ...
               'Temp_n_pcm=0,047619_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=1000_N3=50.mat'};

%
fig1 = figure(20); clf; ax1=gca; set(ax1, 'YScale', 'log'); hold on;
fig2 = figure(21); clf; ax2=gca; set(ax2, 'YScale', 'log'); hold on;

           
for i=1:length(T_test_dirs)
    
    T_test_data = load([root_dir, T_test_dirs{i}]);
    T_test = T_test_data.T(:,T_test_data.N1 + T_test_data.N3);
    
    relErr = abs(1 - T_test ./ T_ref);
    

    plot(ax1, T_test_data.T_ref_dsc, relErr, 'Linewidth', 1.3, 'DisplayName', sprintf('N_{Ag}=%s', num2str(T_test_data.N1)));
    
    %T_test_data.N1
    
    spatial_grid = build_grid(T_test_data.N1, T_test_data.N3, T_test_data.L1, T_test_data.L3, T_test_data.n_tr, T_test_data.n_m, T_test_data.t);
    
    plot(ax2, spatial_grid, 'x', 'DisplayName', sprintf('N_{Ag}=%s', num2str(T_test_data.N1)));
    
end

legend(ax1, 'show', 'location', 'southwest')
xlabel(ax1, 'T_{ref} [degC]')
ylabel(ax1, '|Relative error|')
set(ax1,'FontSize',12);
set(fig1, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.43, 0.6]); 
print(fig1, '/home/argo/masterarbeit/simulationen-data/grid_error/mod_N1_relErr', '-dpng', '-r200');


legend(ax2, 'show', 'location', 'northeast')
xlabel('$$\tilde{x}$$, Computation grid', 'Interpreter', 'latex')
ylabel('$$\frac{\partial \chi}{\partial \tilde{x}}(\tilde{x}) = \Delta x$$ [mm]', 'Interpreter', 'latex')
set(ax2,'FontSize',12);
set(fig2, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.43, 0.6]); 
print(fig2, '/home/argo/masterarbeit/simulationen-data/grid_error/mod_N1_gridsize', '-dpng', '-r200');






%%
%  modify N

root_dir = '/home/argo/masterarbeit/simulationen-data/grid_error/';

T_ref_file = [root_dir, 'Temp_N1=20000_N3=50_L1=40.00_L3=0.10_relTol=1e-7_new.mat'];
% T_ref_file = [root_dir, 'Temp_N1=20000_N3=50_L1=40.00_L3=0.10_relTol=1e-8.mat'];

T_ref_data = load(T_ref_file);
T_ref = T_ref_data.T(:,20000);

T_test_dirs = {'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=160_N3=40.mat', ...
               'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=240_N3=60.mat', ...
               'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=400_N3=100.mat', ...
               'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=800_N3=200.mat', ...
               'Temp_n_pcm=0,2_n_tr=0,1_n_m=0,01_t=0,999_L1=40_L3=0,1_N1=1600_N3=400.mat'};

fig1 = figure(20); clf; ax1=gca; set(ax1, 'YScale', 'log'); hold on;
fig2 = figure(21); clf; ax2=gca; set(ax2, 'YScale', 'log'); hold on;

           
for i=1:length(T_test_dirs)
  
    T_test_data = load([root_dir, T_test_dirs{i}]);
    T_test = T_test_data.T(:,T_test_data.N1);
    
    relErr = abs(1 - T_test ./ T_ref);
%     relErr = (1 - T_test ./ T_ref);
    
    plot(ax1, T_test_data.T_ref_dsc, relErr, 'Linewidth', 1.3, 'DisplayName', sprintf('N=%s', num2str(T_test_data.N1 + T_test_data.N3)));
    
    spatial_grid = build_grid(T_test_data.N1, T_test_data.N3, T_test_data.L1, T_test_data.L3, T_test_data.n_tr, T_test_data.n_m, T_test_data.t);
    
    plot(ax2, spatial_grid, 'x', 'DisplayName', sprintf('N=%s', num2str(T_test_data.N1 + T_test_data.N3)));
    
end

legend(ax1, 'show', 'location', 'southwest')
xlabel(ax1, 'T_{ref} [degC]')
ylabel(ax1, '|Relative error|')
set(ax1,'FontSize',12);
set(fig1, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.43, 0.6]); 
print(fig1, '/home/argo/masterarbeit/simulationen-data/grid_error/mod_N_relErr', '-dpng', '-r200');


legend(ax2, 'show', 'location', 'northeast')
xlabel('$$\tilde{x}$$, Computation grid', 'Interpreter', 'latex')
ylabel('$$\frac{\partial \chi}{\partial \tilde{x}}(\tilde{x}) = \Delta x$$ [mm]', 'Interpreter', 'latex')
set(ax2,'FontSize',12);
set(fig2, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.43, 0.6]); 
print(fig2, '/home/argo/masterarbeit/simulationen-data/grid_error/mod_N_gridsize', '-dpng', '-r200');



%% Relative Error for integration tolerance 1e-7 <-> 1e-8

% Used c_p parametrization with Fraser-Suzuki-Peak, analog in grid error plots
%h  =  10.0; r = 2.0; wr = 5.0; sr = 0.3; z = 130.0; m = 0.3; b = 2.0;

dsc_filename = 'ExpDat_16-407-3_mitKorr_5Kmin_H.csv';
dsc = DSC204_readFile(dsc_filename);
index_T_dsc = [find(dsc.data(:,1) > 29, 1, 'first'), ...
               find(dsc.data(:,1) < 157.9, 1, 'last')];
T_ref_dsc = dsc.data(index_T_dsc(1):index_T_dsc(2),1);

root_dir = '/home/argo/masterarbeit/simulationen-data/grid_error/';

T_ref_file7 = [root_dir, 'Temp_equidistant_grid_N1=20000_N3=50_L1=40,00_L3=0,10_1e-7.mat'];  % 1e-7
% T_ref_file7 = [root_dir, 'Temp_N1=20000_N3=50_L1=40.00_L3=0.10_relTol=1e-7_new.mat'];
T_ref_data7 = load(T_ref_file7);
T_ref7 = T_ref_data7.T(:,20000);

T_ref_file8 = [root_dir, 'Temp_N1=20000_N3=50_L1=40.00_L3=0.10_relTol=1e-8.mat'];    % 1e-8
T_ref_data8 = load(T_ref_file8);
T_ref8 = T_ref_data8.T(:,20000);

relErr = abs(1 - T_ref7 ./ T_ref8);

fig1 = figure(23); clf; ax1 = gca; set(ax1, 'YScale', 'log'); hold on

plot(ax1, T_ref_dsc, relErr, 'Linewidth', 1.3);
xlabel(ax1, 'T_{ref} [degC]')
ylabel(ax1, 'Absolute relative error')
set(ax1,'FontSize',12);
print(fig1, '/home/argo/masterarbeit/thesis/images/integration_tolerance_relErr', '-dpng', '-r200');




%% Get average value of measured heat rate

% dsc_filename = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';
% dsc = DSC204_readFile(dsc_filename);
% 
% idx_start = find(dsc.data(:,2) > 45, 1, 'first');
% heat_rate_meas = gradient(dsc.data(idx_start:end,1), dsc.data(idx_start:end,2)); % [K/min]
% %plot(dsc.data(:,2), dsc.data(:,1));
% plot(dsc.data(idx_start:end,2), heat_rate_meas)
% 
% mean(heat_rate_meas)

dsc_filename = 'ExpDat_16-407-3_mitKorr_2,5Kmin_H.csv';
dsc = DSC204_readFile(dsc_filename);

idx_start = 1;
idx_start = find(dsc.data(:,2) > 150, 1, 'first');
idx_end = find(dsc.data(:,2) < 201, 1, 'last');


heat_rate_meas = gradient(dsc.data(idx_start:idx_end,1), dsc.data(idx_start:idx_end,2));
heat_rate_meas_dev = 100*(heat_rate_meas ./ dsc.Tinfo.Tstep - 1);

fig = figure(50); clf; hold on
plot(dsc.data(idx_start:idx_end,1), heat_rate_meas, 'color', 'blue', ...
    'DisplayName', 'Heat rate at reference crucible');
plot([dsc.data(idx_start,1), dsc.data(end,1)], mean(heat_rate_meas)*ones(1,2), ...
    'Color', 'red', 'Linewidth', 1.3, 'Linestyle', '-', 'DisplayName', 'Mean value');
yticks([2.48, 2.485, 2.49, 2.495, 2.5, 2.505, 2.51]);
yticklabels({'2.48 (-0.8%)', '2.485 (-0.6%)', '2.49 (-0.4%)', '2.495 (-0.2%)', ...
             '2.5 (0%)', '2.505 (+0.2%)', '2.51 (+0.4%)'})


% plot(dsc.data(idx_start:idx_end,2), dsc.data(idx_start:idx_end,3))
% plot(dsc.data(idx_start:idx_end,2), dsc.data(idx_start:idx_end,1))

ylim([2.48 2.511])
xlabel('T_{ref} [degC]');
ylabel('Heat rate (deviation) [K/min]')
legend('show', 'location', 'northwest')
set(gca,'FontSize',12);


print(fig, '/home/argo/masterarbeit/thesis/images/heat_rate_ref_crucible', '-dpng', '-r200');


mean(heat_rate_meas)
var(heat_rate_meas)


%% Berechnung Konfidenzintervalle der geschaetzten Parameter

% fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-10_14:33:40_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS/';
fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-08_22:22:31_407_L1=40_L3=0,1_N1=300_N3=50_5Gaussians/';


file_list = dir(fit_dir);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];

for j=1:length(nameSubDirs)
    
    fit_data = load([fit_dir, nameSubDirs{j}, '/fit_data.mat']);
    p_optim_scaled = fit_data.optimization.p_optim_end;
    param_type = fit_data.optimization.c_p_param_type;
    
    % transform to column vector
    p_optim_scaled = reshape(p_optim_scaled, length(p_optim_scaled),1);
    
    p_optim_unscaled = reverse_scale_params(p_optim_scaled, fit_data.optimization.c_p_param_type);    
    p_optim_unscaled = p_optim_unscaled(fit_data.optimization.p_optim_estimable);
    
    J1_scaled = fit_data.optimization.dqdp_end(:,fit_data.optimization.p_optim_estimable);
    J1_unscaled = reverse_scale_J1(J1_scaled, fit_data.optimization.c_p_param_type, ...
                                   fit_data.optimization.p_optim_estimable);
    
    F1 = fit_data.optimization.residuum_end;
    
    confidence_lvl = 0.05;
    [theta, C] = compute_confidence_interval(F1,J1_unscaled,confidence_lvl);
    p_optim_variance = diag(C);
    np = length(theta); 
    
    % Output
    fprintf('Heat rate: %1.2f\n', fit_data.measurement.dsc_data.Tinfo.Tstep);
    
    
%     fprintf('Scaled params:   ');
%     for i=1:np
%         fprintf('%1.3e\t', p_optim_scaled(i));
%     end
    
    if (strcmp(param_type, 'gauss_linear_comb')) 
        fprintf('\nUnscaled params:\n');
        for i=1:5
            fprintf('Gauss%d: ', i);
            fprintf('%+1.3e +/- %1.3e      %+1.3e +/- %1.3e   %+1.3e +/- %1.3e\n', ...
                p_optim_unscaled(3*i-2), theta(3*i-2), ...
                p_optim_unscaled(3*i-1), theta(3*i-1), ...
                p_optim_unscaled(3*i-0), theta(3*i-0))
        end
        fprintf('Linear: %+1.3e +/- %1.3e\n', ...
            p_optim_unscaled(end-1), theta(end-1));
        fprintf('Const:  %+1.3e +/- %1.3e\n', ...
            p_optim_unscaled(end), theta(end));
            
    elseif (strcmp(param_type, 'fraser_suzuki'))
        fprintf('\nUnscaled params: ');
        for i=1:np
            fprintf('  %+1.4e\t', p_optim_unscaled(i));
        end
        
        fprintf('\nCuboid:          ');
        for i=1:np
            fprintf('+/-%1.4e\t', theta(i));
        end

    end
    
%     fprintf('\nVariances:       ');
%     for i=1:np
%         fprintf('%1.3e\t', p_optim_variance(i));
%     end
    

    
    fprintf('\n\n');
    
        
end

%% Berechnung Faktor b, beta der Messdaten-Kovarianzmatrix

fit_data = load('fit_data.mat');
dsc_data = fit_data.measurement.dsc_data;

idx_start = find(dsc_data.data(:,1) > 140, 1, 'first');

q_meas = dsc_data.data(idx_start:end,3) ./ dsc_data.data(idx_start:end,4);

beta_sq = var(q_meas);

fprintf('beta = %1.2e\n', beta_sq);

q_fig = open('q_pcm_in(T_ref).fig');
children = get(gca, 'Children');
q_res = (children(1).YData)';
close;

b_sq = q_res.' * q_res / (length(q_res) - sum(fit_data.optimization.p_optim_estimable));

fprintf('b    = %1.2e\n', b_sq);


%% Residuum aus heat flux plot nachtraeglich zu fit_data.mat hinzufuegen


fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-03_19:49:07_407_L1=40_L3=0.1_N1=500_N3=50/';

file_list = dir(fit_dir);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];

for j=1:length(nameSubDirs)
    
    
    open([fit_dir, nameSubDirs{j}, '/q_pcm_in(T_ref).fig']);
    children = get(gca, 'Children');
    q_res = (children(1).YData)';
    close;
    
    fit_data = load([fit_dir, nameSubDirs{j}, '/fit_data.mat']);
    fit_data.optimization.residuum_end = q_res;
    
    path_fit_data = [fit_dir, nameSubDirs{j}, '/fit_data.mat'];
    save(path_fit_data, '-struct', 'fit_data');

    
end


%% Fuer einen Fit fuer alle Heizraten einen Doppelplot mit heat flux und 
%  c_p erstellen fuer Latex


% fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-08_22:22:31_407_L1=40_L3=0,1_N1=300_N3=50_5Gaussians/';
fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-09_18:33:20_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS/';


file_list = dir(fit_dir);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];


fig = figure(1); clf;
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);

for j=1:length(nameSubDirs)

    clf;
    if (j < 100)
        ax1 = subplot(2,1,1);
        ax2 = subplot(2,1,2);  
        set(fig, 'Units', 'normalized', 'OuterPosition', [0., 0.1, 0.45, 1.]); 
    else
        ax1 = subplot(1,2,1);
        ax2 = subplot(1,2,2);
        set(fig, 'Units', 'normalized', 'OuterPosition', [0., 0., 1.2, 1.]); 
    end
        
    fit_data = load([fit_dir, nameSubDirs{j}, '/fit_data.mat']);
    p_optim_all = fit_data.optimization.p_optim_end;
    
    % c_p(T) plot
    T_plot = 30:0.01:160;
    switch fit_data.optimization.c_p_param_type
        case 'old_atan_formula'
            c_p_plot = c_p_formula(T_plot, p_optim_all(1:6));
        case 'fraser_suzuki'
            c_p_plot = c_p_fs(T_plot, p_optim_all);
        case 'gauss_linear_comb'
            c_p_plot = c_p_gauss_linear_comb(T_plot, p_optim_all);
    end    
    
    plot(ax1, T_plot, c_p_plot, 'DisplayName', 'c_p(T)', 'Linewidth', 2.)
    legend(ax1, 'show', 'location', 'northwest');
    xlabel(ax1, 'T [degC]');
    ylabel(ax1, 'c_p [mJ/(mg*K)]');
    set(ax1,'FontSize',15, 'FontWeight', 'bold');
    xlim(ax1, [30, 160]);
    
    
    % heat flux plot
    index_T_dsc = fit_data.measurement.index_T_dsc;
    meas_data = fit_data.measurement.dsc_data.data;
    T_ref_dsc = meas_data(index_T_dsc(1):index_T_dsc(2), 1);
    
    q_res = fit_data.optimization.residuum_end;
    q_meas = fit_data.measurement.dsc_data.mass * ...
             meas_data(index_T_dsc(1):index_T_dsc(2),3) ./ ...
             meas_data(index_T_dsc(1):index_T_dsc(2),4);
    q_sim = q_res + q_meas;
    
    plot(ax2, T_ref_dsc, q_sim, 'DisplayName', 'Simulation', 'Linewidth', 2.); hold on
    plot(ax2, T_ref_dsc, q_meas, 'DisplayName', 'Measurement', ...
         'Linestyle', '--', 'Linewidth', 2.);
    plot(ax2, T_ref_dsc, q_res, 'DisplayName', 'Residuum', 'Linewidth', 2.);
    legend(ax2, 'show', 'location', 'northwest');
    xlabel(ax2, 'T_{ref} [degC]');
    ylabel(ax2, '\Phi_q^{PCM,in} [mW]');
    set(ax2,'FontSize',15, 'FontWeight', 'bold');
    xlim(ax2, [T_ref_dsc(1), T_ref_dsc(end)]);

    print(fig, [fit_dir, nameSubDirs{j}, '/combined_img'], '-dpng', '-r300');
    
end
% close;


%% Fuer alle Heizraten c_p(T) und heat flux seperat plotten als .png
%  c_p erstellen fuer Latex


fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-08_22:22:31_407_L1=40_L3=0,1_N1=300_N3=50_5Gaussians/';
% fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-09_18:33:20_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS/';


file_list = dir(fit_dir);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];


fig1 = figure(1); clf; ax1 = gca; hold on;
set(fig1, 'Units', 'normalized', 'OuterPosition', [0., 0., 0.8, 1.]); 
fig2 = figure(2); clf; ax2 = gca; hold on;
set(fig2, 'Units', 'normalized', 'OuterPosition', [0., 0., 0.8, 1.]); 

for j=1:length(nameSubDirs)

    cla(ax1);
    cla(ax2);
    
    fit_data = load([fit_dir, nameSubDirs{j}, '/fit_data.mat']);
    p_optim_all = fit_data.optimization.p_optim_end;
    
    % c_p(T) plot
    T_plot = 30:0.01:160;
    switch fit_data.optimization.c_p_param_type
        case 'old_atan_formula'
            c_p_plot = c_p_formula(T_plot, p_optim_all(1:6));
        case 'fraser_suzuki'
            c_p_plot = c_p_fs(T_plot, p_optim_all);
        case 'gauss_linear_comb'
            c_p_plot = c_p_gauss_linear_comb(T_plot, p_optim_all);
    end    
    
    plot(ax1, T_plot, c_p_plot, 'DisplayName', 'c_p(T)', 'Linewidth', 2.)
    legend(ax1, 'show', 'location', 'northwest');
    xlabel(ax1, 'T [degC]');
    ylabel(ax1, 'c_p [mJ/(mg*K)]');
    set(ax1,'FontSize',24);
    xlim(ax1, [30, 160]);
    title(ax1, sprintf('Heat rate: %1.2f K/min', fit_data.simulation.heat_rate))

    
    % heat flux plot
    index_T_dsc = fit_data.measurement.index_T_dsc;
    meas_data = fit_data.measurement.dsc_data.data;
    T_ref_dsc = meas_data(index_T_dsc(1):index_T_dsc(2), 1);
    
    q_res = fit_data.optimization.residuum_end;
    q_meas = fit_data.measurement.dsc_data.mass * ...
             meas_data(index_T_dsc(1):index_T_dsc(2),3) ./ ...
             meas_data(index_T_dsc(1):index_T_dsc(2),4);
    q_sim = q_res + q_meas;
    
    plot(ax2, T_ref_dsc, q_sim, 'DisplayName', 'Simulation', 'Linewidth', 2.); hold on
    plot(ax2, T_ref_dsc, q_meas, 'DisplayName', 'Measurement', ...
         'Linestyle', '--', 'Linewidth', 2.);
    plot(ax2, T_ref_dsc, q_res, 'DisplayName', 'Residuum', 'Linewidth', 2.);
    legend(ax2, 'show', 'location', 'northwest');
    xlabel(ax2, 'T_{ref} [degC]');
    ylabel(ax2, '\Phi_q^{PCM,in} [mW]');
    set(ax2,'FontSize',24);
    xlim(ax2, [T_ref_dsc(1), T_ref_dsc(end)]);
    title(ax2, sprintf('Heat rate: %1.2f K/min', fit_data.simulation.heat_rate))
    
    print(fig1, [fit_dir, nameSubDirs{j}, '/c_p'], '-dpng', '-r200');
    print(fig2, [fit_dir, nameSubDirs{j}, '/heat_flux'], '-dpng', '-r200');
    
    
end

close(1);
close(2);




%% Compute average value of resulting c_p(T) distribution


fit_dir = '/home/argo/masterarbeit/fits_data/2017-12-08_22:22:31_407_L1=40_L3=0,1_N1=300_N3=50_5Gaussians/';

file_list = dir(fit_dir);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];

for j=1:length(nameSubDirs)

    
    fit_data = load([fit_dir, nameSubDirs{j}, '/fit_data.mat']);
    p_optim_all_scaled = fit_data.optimization.p_optim_end;
    p_optim_all_unscaled = reverse_scale(p_optim_all_scaled, fit_data.optimization.c_p_param_type);
    
    heat_rate = fit_data.simulation.heat_rate;
    
    % TODO: ueber ganze Verteilung ohne linearen und constanten part
    % integrieren, dann normalisieren und dann davon den Erwartungswert
    % berechnen.
    
    if (strcmp(fit_data.optimization.c_p_param_type, 'gauss_linear_comb'))
        error('Not implemented yet!')
        %average_value = sum(p_optim_all_unscaled(1:3:15) .* p_optim_all_unscaled(3:3:15)) / sum(p_optim_all_unscaled(1:3:15));
    elseif (strcmp(fit_data.optimization.c_p_param_type, 'fraser_suzuki'))
        error('Not implemented yet!')
    end 
        
    fprintf('Heat rate: %1.2f: %1.3f\n', heat_rate, average_value);
    
    
    
    
end


%% Plot optimization progress

fit_dir = '2017-12-10_14:33:40_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS';
fit_path = strcat('/home/argo/masterarbeit/fits_data/', fit_dir, '/');

file_list = dir(fit_path);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];


fig = figure(1); clf;
set(fig, 'Units', 'normalized', 'OuterPosition', [0., 0., 0.5, 1.]); 

ax1 = gca; set(ax1, 'YScale', 'log'); hold on


for i=1:length(nameSubDirs)
        
    cla;
    disp(nameSubDirs{i})
    
    filepath = strcat(fit_path, nameSubDirs{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    num_iterations = length(fit_data.optimization.progress_dx_norm);
    
    plot(ax1, 0:num_iterations, fit_data.optimization.progress_F1_norm, ...
        '-s', 'DisplayName', '||F_1^{(k)}||_2', 'Linewidth', 2.);
    plot(ax1, 0:num_iterations-1, fit_data.optimization.progress_dx_norm, ...
        '-s', 'DisplayName', '||\Deltax^{(k)}||_2', 'Linewidth', 2.);
    plot(ax1, 0:num_iterations-1, fit_data.optimization.progress_t_k, ...
        '-s', 'DisplayName', 't^{(k)}', 'Linewidth', 2.);
    plot(ax1, 0:num_iterations, fit_data.optimization.progress_NOC1, ...
        '-s', 'DisplayName', '||\nablaL^{(k)}||_2', 'Linewidth', 2.);
    
    
    xlim(ax1, [0, num_iterations]);
    xlabel('#Iteration');
    legend(ax1, 'show', 'location', 'southwest');
    set(ax1,'FontSize',16, 'FontWeight', 'bold');
        
    print(fig, [fit_path, nameSubDirs{i}, '/optimization_progress'], '-dpng', '-r300');

    
end


%% Compute Melting enthalpy from fit_data

fit_dir = '2017-12-10_14:33:40_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS';
% fit_dir = '2017-12-08_22:22:31_407_L1=40_L3=0,1_N1=300_N3=50_5Gaussians';

fit_path = strcat('/home/argo/masterarbeit/fits_data/', fit_dir, '/');

file_list = dir(fit_path);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];

for i=1:length(nameSubDirs)
        
    filepath = strcat(fit_path, nameSubDirs{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    p_optim_end = fit_data.optimization.p_optim_end;
    
    p_optim_end(end-1:end) = 0.;  
    % turn off linear and constant part, s.t. we just have the phase transition part
    
    switch fit_data.optimization.c_p_param_type
        case 'gauss_linear_comb'
            c_p_fct_handle = @(T) c_p_gauss_linear_comb(T, p_optim_end);
        case 'fraser_suzuki'
            c_p_fct_handle = @(T) c_p_fs(T, p_optim_end);
    end
    
    dH = integral(c_p_fct_handle, 30,200);
    
    fprintf('Heat rate: %1.2g\t dH = %1.2f\n',fit_data.simulation.heat_rate, dH);
    
    
    
end

fprintf('\n');



%% Plot 1st order optimality von lsqnonlin

fit_dir = '2017-12-11_16:26:52_407_L1=40_L3=0.1_N1=300_N3=50';
fit_path = strcat('/home/argo/masterarbeit/fits_data/', fit_dir, '/');

file_list = dir(fit_path);

isub = [file_list(:).isdir]; %# returns logical vector
nameSubDirs = {file_list(isub).name}';
nameSubDirs(ismember(nameSubDirs,{'.','..'})) = [];

fig1 = figure();
ax1 = gca; set(ax1, 'YScale', 'log'); hold on;

for i=1:length(nameSubDirs)
        
    cla;
    
    filepath = strcat(fit_path, nameSubDirs{i}, '/fit_data.mat');
    fit_data = load(filepath);
    
    plot(ax1, fit_data.optimization.firstorderopt, ...
        'DisplayName', num2str(fit_data.simulation.heat_rate));
    
    legend(ax1, 'show', 'location', 'northeast')
    
    
    print(fig1, [fit_path, nameSubDirs{i}, '/optimization_progress_NOC1'], '-dpng', '-r200');
    
end


%% Gauss plot mit allen Gaussen einzeln

fit_data = load('/home/argo/masterarbeit/fits_data/2017-12-08_22:22:31_407_L1=40_L3=0,1_N1=300_N3=50_5Gaussians/2017-12-08_23:53:54_407_0,3Kmin_L1=40_L3=0,1/fit_data.mat');
scaled_params = fit_data.optimization.p_optim_end.';
params = reverse_scale_params(scaled_params, fit_data.optimization.c_p_param_type);

fig1 = figure(1); clf;
ax1 = gca; hold on

T_domain = 30:0.01:160;

m = params(end-1);
b = params(end);

for i=1:5
    
    A_i = params(3*i-2);
    var_i = params(3*i-1);
    offset_i = params(3*i-0);
    
    gauss_eval = A_i * exp(-(T_domain-offset_i).^2 ./ var_i) + m*T_domain + b;
    
    plot(ax1, T_domain, gauss_eval, 'DisplayName', sprintf('Gauss%d', i), ...
        'Linestyle', '-.', 'Linewidth', 1.3);
    
end

plot(ax1, T_domain, c_p_gauss_linear_comb(T_domain, scaled_params), ...
    'DisplayName', 'Total c_p', 'Linewidth', 1.3, 'color', 'black');

set(gca,'FontSize',12)
xlabel(ax1, 'T [degC]')
ylabel(ax1, 'c_p [mJ/(mg*K)]');
legend(ax1, 'show', 'location', 'northwest')
xlim(ax1, [100 140]);

save_path = '/home/argo/masterarbeit/thesis/images/c_p_example'
print(fig1, save_path, '-dpng', '-r200');
    

%% Comparison FwdSens with/without error control

root_dir = '/home/argo/masterarbeit/';

fwdSens_errTol_inactive_path = [root_dir, 'dqdp_gauss_linear_comb_fwdSensTol=inactive.mat'];
fwdSens_errTol_active_path = [root_dir, 'dqdp_gauss_linear_comb_fwdSensTol=1e-7.mat'];

fwdSens_errTol_inactive_data = load(fwdSens_errTol_inactive_path);
fwdSens_errTol_inactive = fwdSens_errTol_inactive_data.Jac;
% fwdSens_errTol_inactive(fwdSens_errTol_inactive > 1e-8) = nan;

fwdSens_errTol_active_data = load(fwdSens_errTol_active_path);
fwdSens_errTol_active = fwdSens_errTol_active_data.Jac;
% fwdSens_errTol_active(fwdSens_errTol_active > 1e-8) = nan;

relErr = 1 - fwdSens_errTol_inactive ./ fwdSens_errTol_active;

% fwdSens unter gewisser Schranke abschneiden weil sonst numerischer
% relativer Fehler keinen Sinn macht. 1e-2 ist aber eig. zu hoch und selbst
% hier ist das Max. bei 12 ...
relErr(fwdSens_errTol_active < 1e-2) = 0;
relErr(fwdSens_errTol_inactive < 1e-2) = 0;


image(relErr(:,[1:9, 31]), 'CDataMapping', 'scaled');
colorbar;

return

figure(1)
image(fwdSens_errTol_inactive, 'CDataMapping', 'scaled');
colorbar;

figure(2)
image(fwdSens_errTol_active, 'CDataMapping', 'scaled');
colorbar;



