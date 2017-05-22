f=1;
N1 = 1000*f;
L1 = 25.;
N2 = 200;
L2 = 1.;
lambda = 0.96;

T_end = 300;

% PCM side
common_args = {'N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, 'lambda', lambda, 'T_end', T_end};
T_pcm_10 = simulate_1d('heat_rate', 10., common_args{:});
T_pcm_5 = simulate_1d('heat_rate', 5., common_args{:});
T_pcm_1 = simulate_1d('heat_rate', 1., common_args{:});

N2 = 0;
L2 = 0.;

% Reference side
common_args = {'N1', N1, 'L1', L1, 'N2', N2, 'L2', L2, 'lambda', lambda, 'T_end', T_end};
T_ref_10 = simulate_1d('heat_rate', 10., common_args{:});
T_ref_5 = simulate_1d('heat_rate', 5., common_args{:});
T_ref_1 = simulate_1d('heat_rate', 1., common_args{:});

offset = 0.;

% take Delta T between reference and link constantan-pcm and plot against T_ref
fig1 = figure(1);
dT = sum(T_ref_10(:,N1-(f-1):N1), 2) - sum(T_pcm_10(:,N1-(f-1):N1), 2);
plot(T_ref_10(:,N1), dT, 'DisplayName', 'beta=10'); hold on

dT = sum(T_ref_5(:,N1-(f-1):N1), 2) - sum(T_pcm_5(:,N1-(f-1):N1), 2);
plot(T_ref_5(:,N1), dT, 'DisplayName', 'beta=5'); hold on

dT = sum(T_ref_1(:,N1-(f-1):N1), 2) - sum(T_pcm_1(:,N1-(f-1):N1), 2);
plot(T_ref_1(:,N1), dT, 'DisplayName', 'beta=1'); hold on

xlabel('T_{ref}');
ylabel('T_{ref} - T_{pcm-Constantan}');
legend(gca, 'show', 'Location', 'northwest')

%print(fig2, '/home/argo/masterarbeit/simulationen-data/N1_1000_N2_200', '-dpng')


