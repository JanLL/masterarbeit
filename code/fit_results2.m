% In diesem File ist bei den Simulationen rho_pcm = 0.85 fixiert.
% Und es wird Silber statt Constantan benutzt.

%get_fit_results('2017-10-14_19:36:00_407_L1=5_L3=0.5_N1=200_N3=50');

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



print(fig, 'c_p_all', '-dpng', '-r500');
print(fig, 'c_p_all_zoom', '-dpng', '-r500');


print(fig, 'c_p_L1=40_L3=0,01', '-dpng', '-r500');




return







%% Plot heat fluxes for 407 sample
dsc_files = DSC204_readFiles('*407*mitKorr*H.csv');
figure(); hold on
for i=1:7
    
    dsc = dsc_files(i);
    plot(dsc.data(:,1), dsc.data(:,3) / dsc.Tinfo.Tstep, 'DisplayName', num2str(dsc.Tinfo.Tstep));
    
end
legend('show', 'location', 'northwest');



