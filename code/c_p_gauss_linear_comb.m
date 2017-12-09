function [c_p] = c_p_gauss_linear_comb(T, p)
%

scaling_ampl   = [2.3314, 13.3911, -2.6634, 1, 1, 1, 1, 1, 1, 1];
scaling_var    = [473.9, 24.51, 121.8, 30, 30, 30, 30, 30, 30, 30];
scaling_offset = [134.38, 127.9, 145.7, 130, 130, 130, 130, 130, 130, 130];
scaling_linear = 0.01;
scaling_const  = 1.;


c_p = 0;
for i=1:10
    ampl_tmp   = scaling_ampl(i)   * p(3*i - 2);
    var_tmp    = scaling_var(i)    * p(3*i - 1);
    offset_tmp = scaling_offset(i) * p(3*i - 0);
    c_p = c_p + ampl_tmp * exp(-(T-offset_tmp).^2 / var_tmp);
    
end
c_p = c_p + scaling_linear*p(31)*T + scaling_const*p(32);


% TODO: vektorisieren.. momentan problem mit dimensionen wenn T vektor
% Temperatur und parameter in matrizen schreiben und dann ueber richtige
% dimension aufsummieren...
% gauss_index_ampl = 1:3:30;
% gauss_index_sigma = 2:3:30;
% gauss_index_offset = 3:3:30;
% 
% c_p = sum(p(gauss_index_ampl) .* exp(-(T-p(gauss_index_offset)).^2 ./ p(gauss_index_sigma).^2));


end

