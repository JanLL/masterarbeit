function [rho] = PCM_rho( T_K )
% rho [kg/m^3] as function of T [K]
% for Polymer

    % trafo for [K]->[°C]
    T_C = T_K-273.15;

%     T_C = 150;

    c(1) = 167.2182;
    c(2) = 129.9861;
    c(3) = 760.8218;
    c(4) = 0.078916;   
    rho = c(1) ./(1.0+exp(c(4).*(T_C-c(2))  )) +c(3);
    % rho = 935.0; % [kg/m^3]
    
end    

% % data Marlotherm SH datasheet
% T_C = [0 25 50 75 100 120 150 175 ];
% rho = [0.941 0.934 0.925 0.915 0.902 0.885 0.782 0.770 ]*1000;
% plot(T_C, rho,'+')
% 
% x = -10:0.05:10
% y = 1.0./(1+exp(x))
% plot(x,y)
