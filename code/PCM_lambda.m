function [lambda] = PCM_lambda( T_K )
% lambda [W/m/K] as function of T [K]
% for Polymer

%     if (T_K>200+273.15)
%         warning('in < m203_lambda_Polym_T > function only valid up to 200°C')
%     end
    % trafo for [K]->[°C]
    T_C = T_K-273.15;

%     T_C = 150;
    
    % Coefficients:
    c(1) = 0.41857;
    c(2) = 96.162 * 2.762;
    c(3) = 0.15406 * 3.514;
    c(4) = 0.036677;
    lambda = c(1)./(1.0+exp(  c(4).*(T_C-c(2))  )) +c(3);
%     lambda = (c(1)./(1.0+exp(  c(4).*(T_C-c(2))  )) +c(3))*0.0025;

end    

% %   // thermal conductivity of the PCM
% %   parameter Real lambda_pcm=0.3;
% T_C = [0 25 50 75 100 120 150 175 ];
% lambda = [0.57 0.54 0.50 0.43 0.38 0.26 0.19 0.19];
% plot(T_C, lambda,'r+')
% T_C2 = 5:5:200;
% lambda2 = p1*T_C2.^3 + p2*T_C2.^2 + p3*T_C2 + p4;
% hold on
% plot(T_C2, lambda2,'b-')