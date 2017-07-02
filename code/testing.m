datetime_cell = num2cell(int32(clock));
solverName = 'lsqnonlin';
heat_rate = 0.3;
dsc_fileSpec = 'ExpDat_16-407-3_mitKorr_10Kmin_H.csv';


datetime_str = sprintf('%04i-%02i-%02i_%02i:%02i:%02i', datetime_cell{:});

heat_rate_str = num2str(heat_rate);
heat_rate_str = strrep(heat_rate_str, '.', ',');

mass_code_str = dsc_fileSpec(11:13);



final_str = strcat(datetime_str, '_', mass_code_str, '_', heat_rate_str, 'Kmin_', solverName ,'.mat')
