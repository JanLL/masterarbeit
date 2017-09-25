function [stop] = disp_aux(x, optimValues, state)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Save optimization variables for all iteration steps
persistent x_process;
global p_optim_process;

if optimValues.iteration == size(x_process, 1)
    x_process = [x_process; x];
end

if strcmp(state, 'done')
    p_optim_process = x_process;
end


stop = false;

for j=0:ceil(length(x)/10)-1
    fprintf('param number:\t');
    fprintf(repmat('%d\t',1,10),(1:10)+j*10 );
    fprintf('\n');
    
    fprintf('param value:\t');
    index = (1:10)+j*10;
    index = index(index <= length(x));
    fprintf(repmat('%1.3g\t',1,10), x(index));
    fprintf('\n');
end

fprintf('\n');

end

