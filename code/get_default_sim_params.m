function [sim_params_pcm, sim_params_ref] = get_default_sim_params(varargin)
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
%   sim_params --> struct with (default) simulation parameters as input for
%                  simulate_1d.

sim_params_pcm = struct();

% check for input arguments and update variables where necessary
if hasOption(varargin, 'L1'), sim_params_pcm.L1 = getOption(varargin, 'L1'); 
else sim_params_pcm.L1 = 25.; end
if hasOption(varargin, 'L2'), sim_params_pcm.L2 = getOption(varargin, 'L2'); 
else sim_params_pcm.L2 = 0.; end
if hasOption(varargin, 'L3'), sim_params_pcm.L3 = getOption(varargin, 'L3'); 
else sim_params_pcm.L3 = 1.; end

if hasOption(varargin, 'N3'), sim_params_pcm.N3 = getOption(varargin, 'N3');
else sim_params_pcm.N3 = 50; end
if sim_params_pcm.N3 <= 0, error('N3 must be greater than 0!'); end

% compute N1, N2 s.t. dx is equal everywhere as default
if hasOption(varargin, 'N1'), sim_params_pcm.N1 = getOption(varargin, 'N1'); 
else sim_params_pcm.N1 = sim_params_pcm.L1 / sim_params_pcm.L3 * sim_params_pcm.N3; end

if hasOption(varargin, 'N2'), sim_params_pcm.N2 = getOption(varargin, 'N2'); 
else sim_params_pcm.N2 = sim_params_pcm.L2 / sim_params_pcm.L3 * sim_params_pcm.N3; end

if hasOption(varargin, 'heat_rate'), sim_params_pcm.heat_rate = getOption(varargin, 'heat_rate'); 
else sim_params_pcm.heat_rate = 10; end
if hasOption(varargin, 'T_0'), sim_params_pcm.T_0 = getOption(varargin, 'T_0'); 
else sim_params_pcm.T_0 = 10.; end
if hasOption(varargin, 'T_end'), sim_params_pcm.T_end = getOption(varargin, 'T_end'); 
else sim_params_pcm.T_end = 300.; end

sim_params_ref = deal(sim_params_pcm);
sim_params_ref.N3 = 0;
sim_params_ref.L3 = 0.;


end

