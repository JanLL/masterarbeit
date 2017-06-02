function [sim_params] = get_param_sim(varargin)
%
%
% INPUT (variable):          
%        N1 --> number of spatial discretization lattice points (Constantan)
%        N2 --> number of spatial discretization lattice points (Crucible)
%        N3 --> number of spatial discretization lattice points (PCM)
%        L1 --> length [mm] of Constantan part.
%        L2 --> length [mm] of crucible part.
%        L3 --> length [mm] of PCM part.
%       T_0 --> initial temperature [degree celsius] everywhere.
%     T_end --> integrate until T_oven as reached T_end [degree celsius].
% heat_rate --> rate [K/min] with which oven temperature increases.
%
% OUTPUT: 
%   sim_params --> 1x2 struct with (default) simulation parameters as input 
%                  for simulate_1d. 
%                  sim_params(1): simulation with PCM
%                  sim_params(2): reference simulation without PCM

sim_params = struct();

% check for input arguments and update variables where necessary
if hasOption(varargin, 'L1'), sim_params.L1 = getOption(varargin, 'L1'); 
else sim_params.L1 = 25.; end
if hasOption(varargin, 'L2'), sim_params.L2 = getOption(varargin, 'L2'); 
else sim_params.L2 = 0.; end
if hasOption(varargin, 'L3'), sim_params.L3 = getOption(varargin, 'L3'); 
else sim_params.L3 = 1.; end

if hasOption(varargin, 'N3'), sim_params.N3 = getOption(varargin, 'N3');
else sim_params.N3 = 50; end
if sim_params.N3 <= 0, error('N3 must be greater than 0!'); end

% compute N1, N2 s.t. dx is equal everywhere as default
if hasOption(varargin, 'N1'), sim_params.N1 = getOption(varargin, 'N1'); 
else sim_params.N1 = sim_params.L1 / sim_params.L3 * sim_params.N3; end

if hasOption(varargin, 'N2'), sim_params.N2 = getOption(varargin, 'N2'); 
else sim_params.N2 = sim_params.L2 / sim_params.L3 * sim_params.N3; end

if hasOption(varargin, 'heat_rate'), sim_params.heat_rate = getOption(varargin, 'heat_rate'); 
else sim_params.heat_rate = 10; end
if hasOption(varargin, 'T_0'), sim_params.T_0 = getOption(varargin, 'T_0'); 
else sim_params.T_0 = 10.; end
if hasOption(varargin, 'T_end'), sim_params.T_end = getOption(varargin, 'T_end'); 
else sim_params.T_end = 300.; end

sim_params(2) = deal(sim_params(1));
sim_params(2).N3 = 0;
sim_params(2).L3 = 0.;

end

