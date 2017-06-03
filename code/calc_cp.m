function [c_p_pcm] = calc_cp()

dsc = DSC204_readFile('ExpDat_16-407-3_mitKorr_10Kmin_H.csv');
dsc_sap = DSC204_readFile('Waermekapazitaet_Saphirmessung/Sap-Kurve_10Kmin_H_Segment_7.csv');

% masses actually not needed because voltage already normalized by mass:
% U_korr_pcm = (U_pcm - U_0)/m_pcm}
%m_pcm = dsc.mass;
%m_sap = dsc_sap.mass;

% we need to interpolate the data on the same temp-grid because
% dsc.data(:,1) != dsc_sap.data(:,1)
T_grid = dsc.data(:,1);
U_korr_pcm = dsc.data(:,3);
U_korr_sap = interp1(dsc_sap.data(:,1), dsc_sap.data(:,3), T_grid);

c_p_sap = DSC204_cp_saphire_DIN11357(T_grid);

c_p_pcm = c_p_sap .* U_korr_pcm./U_korr_sap; 


end









