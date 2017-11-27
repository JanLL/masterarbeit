function [stop] = disp_aux(x, optimValues, state)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

stop = false;


% Save optimization variables for all iteration steps
persistent x_process;
global p_optim_process;

p_optim_process = false ;



if (optimValues.iteration == size(x_process, 1) && strcmp(state, 'iter'))
    try
        x_process = [x_process; x];
    catch
        size(x_process)
        size(x)
        fprintf('Error with x_process occured!\n');
    end
end

if strcmp(state, 'done')
    p_optim_process = x_process;
end


% Print optimization variables on the command window
if strcmp(state, 'iter')
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
end

fprintf('\n');

end

