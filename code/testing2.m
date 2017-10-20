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

c_p_fig = open('c_p(T).fig');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
children = get(gca, 'Children');
children.LineWidth = 2;
set(gca,'FontSize',24)
set(gca,'xlim', [20 max(children.XData)]);
ylabel('c_p [mJ/(mg*K)]');
title('');

print(c_p_fig, 'c_p(T)', '-dpng', '-r200');
close;

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

print(c_p_fig, 'q_pcm_in(T_ref)', '-dpng', '-r200');
title(title_str);
close;


return


% Read 

dsc_list = DSC204_readFiles('DSC204_F1_Phoenix_Messungen/Messungen/Messungen/ExpDat_16-407-3_mitKorr_*Kmin_H.csv');

figure(5); hold on
for i=1:length(dsc_list)

    dsc = dsc_list(i);
    plot(dsc.data(:,1), dsc.data(:,3) ./ dsc.Tinfo.Tstep, 'DisplayName', num2str(dsc.Tinfo.Tstep));
    
    
end


return






% Convergence Rate Plots
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











% testing of frazer-suzuki parametrization

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
