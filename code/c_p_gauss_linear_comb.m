function [c_p] = c_p_gauss_linear_comb(T, p)
%


c_p_1 = p(1) * exp(-(T-p(3)).^2 / p(2)^2);
c_p_2 = p(4) * exp(-(T-p(6)).^2 / p(5)^2);
c_p_3 = p(7) * exp(-(T-p(9)).^2 / p(8)^2);
c_p_4 = p(10) * exp(-(T-p(12)).^2 / p(11)^2);
c_p_5 = p(13) * exp(-(T-p(15)).^2 / p(14)^2);

c_p = c_p_1 + c_p_2 + c_p_3 + c_p_4 + c_p_5 + p(16); 

end

