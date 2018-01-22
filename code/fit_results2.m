% In diesem File ist bei den Simulationen rho_pcm = 0.85 fixiert.
% Und es wird Silber statt Constantan benutzt.

[fig, enthalpies] = get_fit_results('2017-10-15_21:08:00_407_L1=5_L3=0.1_N1=200_N3=50');


[fig, enthalpies] = get_fit_results('2017-10-17_21:19:17_407_L1=50_L3=0.1_N1=200_N3=50');
[fig, enthalpies] = get_fit_results('2017-10-15_23:58:00_407_L1=60_L3=0.1_N1=200_N3=50');
[fig, enthalpies] = get_fit_results('2017-10-16_04:45:00_407_L1=40_L3=0.1_N1=200_N3=50');


[fig, enthalpies] = get_fit_results('2017-10-16_04:45:00_407_L1=40_L3=0.1_N1=200_N3=50');
[fig, enthalpies] = get_fit_results('2017-10-16_11:18:10_407_L1=40_L3=0.13_N1=200_N3=50');
[fig, enthalpies] = get_fit_results('2017-10-16_10:08:10_407_L1=40_L3=0.07_N1=200_N3=50');

[fig, enthalpies] = get_fit_results('2017-10-18_00:17:12_407_L1=30_L3=0.1_N1=500_N3=50');

[fig, enthalpies] = get_fit_results('2017-10-18_10:59:34_407_L1=80_L3=0.1_N1=200_N3=50');


[fig, enthalpies] = get_fit_results('2017-10-18_12:11:40_407_L1=40_L3=0.01_N1=200_N3=50');
[fig, enthalpies] = get_fit_results('2017-10-18_13:28:22_407_L1=5_L3=0.01_N1=200_N3=50');


[fig, enthalpies] = get_fit_results('2017-10-18_14:41:19_407_L1=40_L3=0.2_N1=200_N3=50');


[fig, enthalpies] = get_fit_results('2017-10-18_16:40:21_407_L1=40_L3=0.1_N1=500_N3=50');

[fig, enthalpies] = get_fit_results('2017-10-23_21:04:27_407_L1=40_L3=0.1_N1=500_N3=50');

[fig, enthalpies] = get_fit_results('2017-10-29_20:40:04_407_L1=40_L3=0.1_N1=500_N3=50');

% Fraser-Suzuki Peak
[fig, enthalpies] = get_fit_results('2017-11-08_16:45:36_407_L1=40_L3=0.1_N1=500_N3=50');


% Modified heat rate
[fig, enthalpies] = get_fit_results('2017-11-27_08:22:19_407_L1=40_L3=0.1_N1=500_N3=50'); % FS
[fig, enthalpies] = get_fit_results('2017-11-27_17:36:31_407_L1=40_L3=0.1_N1=500_N3=50'); % Gausse

[fig, enthalpies] = get_fit_results('2017-11-27_22:04:37_407_L1=40_L3=0.1_N1=500_N3=200'); % FS

% Unmodified heat rate again
[fig, enthalpies] = get_fit_results('2017-12-03_19:49:07_407_L1=40_L3=0.1_N1=500_N3=50'); % Gausse


% Own GN solver insteadt lsqnonlin
% Nominal heat rate
[fig, enthalpies] = get_fit_results('2017-12-07_17:07:33_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS');  % nominal heat rate
[fig, enthalpies] = get_fit_results('2017-12-08_15:36:28_407_L1=40_L3=0,1_N1=300_N3=50_GN_Gaussians');  % 3 Gaussians
[fig, enthalpies] = get_fit_results('2017-12-08_19:30:45_407_L1=40_L3=0,1_N1=300_N3=50_GN_Gaussians');  % 5 Gaussians


[fig, enthalpies] = get_fit_results('2017-12-08_22:22:31_407_L1=40_L3=0,1_N1=300_N3=50_5Gaussians');  % 5 Gaussians
[fig, enthalpies] = get_fit_results('2017-12-09_11:50:41_407_L1=40_L3=0.1_N1=300_N3=50_3Gaussians');  % 3 Gaussians


[fig, enthalpies] = get_fit_results('2017-12-09_18:33:20_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS');  % FS
[fig, enthalpies] = get_fit_results('2017-12-09_22:35:53_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS_modHeatRate');  % FS


[fig, enthalpies] = get_fit_results('2017-12-10_11:59:41_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS');  % FS
[fig, enthalpies] = get_fit_results('2017-12-10_14:33:40_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS');  % FS

% test with n_tr = 0.99 instead 0.999
[fig, enthalpies] = get_fit_results('2017-12-11_18:35:22_407_L1=40_L3=0.1_N1=300_N3=50');  % FS

% fwdSensTol decreased
[fig, enthalpies] = get_fit_results('2017-12-14_23:47:55_407_L1=40_L3=0,1_N1=200_N3=50_GN_FS');  % FS

% finally used ones
[fig, enthalpies] = get_fit_results('2017-12-20_14:25:10_407_L1=40_L3=0,1_N1=300_N3=50_GN_FS_used');  % FS
[fig, enthalpies] = get_fit_results('2017-12-19_20:27:59_407_L1=40_L3=0.1_N1=300_N3=50_5Gaussians_used');  % FS




figure(66);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.66 1]);
children = get(gca, 'Children');
num_children = length(children);
for i=1:num_children; children(i).LineWidth = 1.3; end
set(gca,'FontSize',24)
set(gca,'xlim', [105 160]);
set(gca,'ylim', [0 60]);
ylabel('c_p [mJ/(mg*K)]');
title('Linear combination of Gaussians');
% title('Fraser-Suzuki peak')
print(fig, 'c_p_all_Gaussians', '-dpng', '-r200');
% print(fig, 'c_p_all_zoom_Gaussians', '-dpng', '-r200');
% savefig(fig, 'c_p_all.fig');

% set(gca,'xlim', [110 160]);
% %set(gca,'ylim', [0 330]);
% print(fig, 'c_p_all_zoom', '-dpng', '-r200');
close();



% print(fig, 'c_p_L1=40_L3=0,01', '-dpng', '-r200');




return







%% Plot heat fluxes for 407 sample
dsc_files = DSC204_readFiles('*417*mitKorr*H.csv');
figure(); hold on
for i=1:7
    
    dsc = dsc_files(i);
    plot(dsc.data(:,1), dsc.data(:,3) / dsc.Tinfo.Tstep, 'DisplayName', num2str(dsc.Tinfo.Tstep));
    
end
legend('show', 'location', 'northwest');



