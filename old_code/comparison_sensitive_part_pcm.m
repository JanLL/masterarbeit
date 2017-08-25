% path for data extraction functions of measurements
path('/home/argo/masterarbeit/DSC204_F1_Phoenix_Messungen', path);

dsc_list = DSC204_readFiles(['/home/argo/masterarbeit/DSC204_F1_Phoenix_Messungen/' ...
    'Messungen/Messungen/ExpDat_16-407-3_mitKorr_*Kmin_H.csv']);

% skipped 0.3 and 0.6 K/min because no data in sensitive area
for i=3:length(dsc_list)
    dsc = dsc_list(i);
    
    T_domain = dsc.data(:,1);
    U_domain = dsc.data(:,3);
    
    %plot(T_domain, U_domain, 'DisplayName',strcat(num2str(dsc.Tinfo.Tstep), dsc.Tinfo.Tstepunit));
    %hold on
    
    indices_sensitive = find(dsc_list(i).data(:,1)>30 & dsc_list(i).data(:,1)<80);
    T_domain = T_domain(indices_sensitive);
    U_domain = U_domain(indices_sensitive);
    
    plot(T_domain, U_domain ./ dsc.Tinfo.Tstep, ...
        'DisplayName',strcat(num2str(dsc.Tinfo.Tstep), dsc.Tinfo.Tstepunit));
    hold on
    
    % fit linear function to get slope c
    linear_fit = fittype('a*x + b', 'coeff', {'a', 'b'});
    start_param = [1., 0.];
    f = fit(T_domain, U_domain ./ dsc.Tinfo.Tstep, linear_fit, ...
        'Startpoint', start_param);
    
    coeffs_average = coeffvalues(f);
    slope_average = coeffs_average(1);
    
    coeffs_conf = confint(f);
    slope_std_dev = coeffs_conf(2,1) - slope_average;
    
    disp(sprintf('%e +- %e', slope_average, slope_std_dev));
    
    % Mittelwert ist NICHT im 1-sigma bereich der anderen werte jeweils 
    % aber: als Fehler wurde hier nur der Fit-Fehler genommen und keine
    % Messfehler, welche sehr wahrsch. deutlich groesser sind.
end

legend('show', 'Location','northwest');
hold off