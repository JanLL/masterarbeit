startup;

N1 = 50;
L1 = 50.;
N2 = 450;
L2 = 5.;
lambda = 0.96;

T_end = 1000;


T_pcm_10 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                   'lambda', lambda, 'heat_rate', 10., 'T_end', T_end);
T_pcm_5 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                  'lambda', lambda, 'heat_rate', 5., 'T_end', T_end);
T_pcm_1 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                  'lambda', lambda, 'heat_rate', 1., 'T_end', T_end);

N2 = 0;
L2 = 0.;

T_ref_10 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                   'lambda', lambda, 'heat_rate', 10., 'T_end', T_end);
T_ref_5 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                  'lambda', lambda, 'heat_rate', 5., 'T_end', T_end);
T_ref_1 = simulate_1d('N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, ...
                  'lambda', lambda, 'heat_rate', 1., 'T_end', T_end);
     
% take Delta T between reference and end of pcm and plot against T_ref
figure(1)
plot(T_ref_1(:,N1), (T_ref_1(:,N1)-T_pcm_1(:,end))*1., 'DisplayName', 'beta=1'); hold on
plot(T_ref_5(:,N1), (T_ref_5(:,N1)-T_pcm_5(:,end))*1., 'DisplayName', 'beta=5'); hold on
plot(T_ref_10(:,N1), (T_ref_10(:,N1)-T_pcm_10(:,end)), 'DisplayName', 'beta=10'); hold on

axis([70, 200, -inf, inf]);
xlabel('T_{ref}');
ylabel('T_{ref} - T_{pcm,N2}');
legend(gca, 'show', 'Location', 'northeast')

% take Delta T between reference and link constantan-pcm and plot against T_ref
figure(2)
plot(T_ref_1(:,N1), (T_ref_1(:,N1)-T_pcm_1(:,N1))*1., 'DisplayName', 'beta=1'); hold on
plot(T_ref_5(:,N1), (T_ref_5(:,N1)-T_pcm_5(:,N1))*1., 'DisplayName', 'beta=5'); hold on
plot(T_ref_10(:,N1), (T_ref_10(:,N1)-T_pcm_10(:,N1)), 'DisplayName', 'beta=10'); hold on

axis([70, 200, -inf, inf]);
xlabel('T_{ref}');
ylabel('T_{ref} - T_{pcm,N2}');
legend(gca, 'show', 'Location', 'northeast')

% comparison of simulation and measurements
dsc = DSC204_readFile(['/home/argo/masterarbeit/DSC204_F1_Phoenix_Messungen/' ...
    'Messungen/Messungen/ExpDat_16-407-3_mitKorr_10Kmin_H.csv']);

figure(3)
plot(T_ref_10(:,N1), (T_ref_10(:,N1)-T_pcm_10(:,end)) / max(T_ref_10(:,N1)-T_pcm_10(:,end)) * max(dsc.data(:,3)), 'DisplayName', 'simulation, beta=10.'); hold on
plot(dsc.data(:,1), dsc.data(:,3), 'DisplayName', 'measurement, beta=10.'); hold on

axis([70, 200, -inf, inf])
legend(gca, 'show', 'Location', 'northeast');
title('normalized on maximum value');
xlabel('T_{ref}');
ylabel('\Delta T');
legend(gca, 'show', 'Location', 'northeast')


