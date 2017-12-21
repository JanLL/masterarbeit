function c_p_meas = calc_cp(dsc_pcm_src)
% TODO: Description!!

if ischar(dsc_pcm_src)
    dsc_pcm = DSC204_readFile(dsc_pcm_src);
elseif isstruct(dsc_pcm_src)
    dsc_pcm = dsc_pcm_src;
else
    error('dsc_pcm_src must be either string with path or struct of dsc_pcm data!');
end


dsc_sap = DSC204_readFile('Waermekapazitaet_Saphirmessung/Sap-Kurve_10Kmin_H_Segment_7.csv');

heat_rate_pcm = dsc_pcm.Tinfo.Tstep;
heat_rate_sap = 10.;

heat_rate_factor = heat_rate_pcm / heat_rate_sap;
% heat_rate_factor because we just have sapphire measurements of 10K/min
% Assuming linear dependency of dU w.r.t. heat_rate (see 
% comparison_sensitive_part_pcm.m for reasoning) we can generate U_korr_sap
% for other heat_rates.

% masses actually not needed because voltage already normalized by mass:
% U_korr_pcm = (U_pcm - U_0)/m_pcm}
%m_pcm = dsc.mass;
%m_sap = dsc_sap.mass;

% we need to interpolate the data on the same temp-grid because
% dsc.data(:,1) != dsc_sap.data(:,1)
T_ref_meas = dsc_pcm.data(:,1);
U_korr_pcm = dsc_pcm.data(:,3);
U_korr_sap = interp1(dsc_sap.data(:,1), heat_rate_factor*dsc_sap.data(:,3), T_ref_meas);

c_p_sap = DSC204_cp_saphire_DIN11357(T_ref_meas);


c_p_meas = c_p_sap .* U_korr_pcm./U_korr_sap; 

c_p_meas = [T_ref_meas, c_p_meas];

end









