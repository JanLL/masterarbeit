function [stop] = disp_aux(x, optimValues, state)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


stop = false;

%disp(x);

% TODO: verallgemeinern fuer beliege groesse von x ... mit 2xlength(x)
% matrix und strcat(s(:,1)', s(:,2)')
x = 1:0.12783213:10;
fprintf('param number:\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n');
fprintf(strcat('param value:\t', repmat('%1.5g\t',1,10)), x(1:10));


end

