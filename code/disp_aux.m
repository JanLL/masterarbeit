function [stop] = disp_aux(x, optimValues, state)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


stop = false;

%disp(x);

% TODO: verallgemeinern fuer beliege groesse von x ... mit 2xlength(x)
% matrix und strcat(s(:,1)', s(:,2)')
x = 1:0.2273:10;

for j=0:ceil(length(x)/10)-1
    fprintf('param number:\t');
    fprintf(repmat('%d\t',1,10),(1:10)+j*10 );
    fprintf('\n');
    
    fprintf('param value:\t');
    index = (1:10)+j*10;
    index = index(index < length(x));
    fprintf(repmat('%1.5g\t',1,10), x(index));
    fprintf('\n');
end

fprintf('\n');

end

