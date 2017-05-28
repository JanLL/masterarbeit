N1 = 1250;
L1 = 25.;
N2 = 0;
L2 = 0.;
N3 = 50;
L3 = 1.;
lambda = 0.96;

T_end = 300;

% PCM side
common_args = {'N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, 'N3', N3, 'L3', L3, ...
               'lambda', lambda, 'T_end', T_end, 'c_p_shape', 'delta_distr'};
T_pcm_10 = simulate_1d('heat_rate', 10., common_args{:});
T_pcm_5 = simulate_1d('heat_rate', 5., common_args{:});
T_pcm_1 = simulate_1d('heat_rate', 1., common_args{:});
% T_pcm_01 = simulate_1d('heat_rate', 0.01, common_args{:});


N3 = 0;
L3 = 0.;

% Reference side
common_args = {'N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, 'N3', N3, 'L3', L3, ...
               'lambda', lambda, 'T_end', T_end, 'c_p_shape', 'delta_distr'};
T_ref_10 = simulate_1d('heat_rate', 10., common_args{:});
T_ref_5 = simulate_1d('heat_rate', 5., common_args{:});
T_ref_1 = simulate_1d('heat_rate', 1., common_args{:});
% T_ref_01 = simulate_1d('heat_rate', 0.01, common_args{:});


offset = 0.;

% take Delta T between reference and link constantan-pcm and plot against T_ref
fig1 = figure(1);
dT10 = T_ref_10(:,N1) - T_pcm_10(:,N1);
plot(T_ref_10(:,N1), dT10, 'DisplayName', 'beta=10'); hold on

dT5 = T_ref_5(:,N1) - T_pcm_5(:,N1);
plot(T_ref_5(:,N1), 2.*dT5, 'DisplayName', 'beta=5'); hold on

dT1 = T_ref_1(:,N1) - T_pcm_1(:,N1);
plot(T_ref_1(:,N1), 10.*dT1, 'DisplayName', 'beta=1'); hold on

% dT01 = T_ref_01(:,N1) - T_pcm_01(:,N1);
% plot(T_ref_01(:,N1), 100.*dT01, 'DisplayName', 'beta=0.1'); hold on

xlabel('T_{ref}');
ylabel('T_{ref} - T_{pcm-Constantan}');
legend(gca, 'show', 'Location', 'northwest')

%print(fig1, '/home/argo/masterarbeit/simulationen-data/Delta_T_delta_peak', '-dpng')


