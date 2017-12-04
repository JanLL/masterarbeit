function [c_p] = c_p_gauss_linear_comb(T, p)
%



% c_p_1 = p(1) * exp(-(T-p(3)).^2 / p(2)^2);
% c_p_2 = p(4) * exp(-(T-p(6)).^2 / p(5)^2);
% c_p_3 = p(7) * exp(-(T-p(9)).^2 / p(8)^2);
% c_p_4 = p(10) * exp(-(T-p(12)).^2 / p(11)^2);
% c_p_5 = p(13) * exp(-(T-p(15)).^2 / p(14)^2);
% 
% c_p = c_p_1 + c_p_2 + c_p_3 + c_p_4 + c_p_5 + p(16)*T + p(17); 

c_p = 0;
for i=1:3:30
    c_p = c_p + p(i) * exp(-(T-p(i+2)).^2 / p(i+1));
    
end
c_p = c_p + 0.01*p(31)*T + p(32);


% TODO: vektorisieren.. momentan problem mit dimensionen wenn T vektor
% Temperatur und parameter in matrizen schreiben und dann ueber richtige
% dimension aufsummieren...
% gauss_index_ampl = 1:3:30;
% gauss_index_sigma = 2:3:30;
% gauss_index_offset = 3:3:30;
% 
% c_p = sum(p(gauss_index_ampl) .* exp(-(T-p(gauss_index_offset)).^2 ./ p(gauss_index_sigma).^2));


end

