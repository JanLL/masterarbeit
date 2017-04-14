% path for data extraction functions of measurements
path('/home/argo/masterarbeit/DSC204_F1_Phoenix_Messungen', path);

dsc = DSC204_readFile(['/home/argo/masterarbeit/DSC204_F1_Phoenix_Messungen/' ...
    'Messungen/Messungen/ExpDat_16-407-3_mitKorr_10Kmin_H.csv']);
dsc_sap = DSC204_readFile(['/home/argo/masterarbeit/DSC204_F1_Phoenix_Messungen/' ...
    'Waermekapazitaet_Saphirmessung/Sap-Kurve_10Kmin_H_Segment_7.csv']);

m_pcm = dsc.mass;
m_sap = dsc_sap.mass;

% we need to interpolate the data on the same temp-grid because
% dsc.data(:,1) != dsc_sap.data(:,1)
T_grid = dsc.data(:,1);
U_korr_pcm = dsc.data(:,3);
U_korr_sap = interp1(dsc_sap.data(:,1), dsc_sap.data(:,3), T_grid);

c_p_sap = DSC204_cp_saphire_DIN11357(T_grid);

c_p_pcm = c_p_sap .* m_sap./m_pcm .* U_korr_pcm./U_korr_sap;
% !!!! mit Massen multipl. weil V/g Einheit U_korr...

% save data for further processing in python
filename = 'c_p_pcm.mat';
save(filename, 'c_p_pcm', 'T_grid');

% fit c_p on function
fit_fct2 = fittype('1/(exp(d*(x-a))+1) * (b*exp(-c*(x-a)^2)) + e*x + f', ...
    'coeff', {'a', 'b', 'c', 'd', 'e', 'f'});
start_param = [140., 100., 0.01, 0.7, 0.01, 15.];
f = fit(T_grid(100:end-10), c_p_pcm(100:end-10), fit_fct2, ...
    'Startpoint', start_param);

c_p_fit = f(T_grid(100:end-10));

plot(T_grid(100:end-10), c_p_fit, 'color', 'red'); hold on
%plot(T_grid(100:end-10), c_p_pcm(100:end-10), 'color', 'blue');
hold off
