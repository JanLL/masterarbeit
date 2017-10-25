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


