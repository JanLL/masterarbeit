function [dU] = load_dU(fit_data_input)
%
%
% INPUT:
%    fit_data_input --> either struct with saved fit data or string with
%                       path to .mat file with saved fit data.
%
% OUTPUT:
%                dU --> (n x 3) array where n number of T_ref evaluation
%                       points. 
%                       1st column: T_ref
%                       2nd column: dU from measurement
%                       3rd column: dU from optimization


if isstruct(fit_data_input)
    fit_data = fit_data_input;
elseif isstr(fit_data_input)
    fit_data = load(fit_data_input);
    
    


end

