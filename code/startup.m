% set root directory
tmpPATHS.code = pwd; 
tmpPATHS.root = fullfile(pwd, '../');

% path to SOLVIND
%tmpPATHS.solvind = fullfile('/home/argo/SOLVIND_SUITE/Packages/SOLVIND/Debug/matlab');
tmpPATHS.solvind = fullfile('/home/argo/mod_SOLVIND_SUITE/Packages/SOLVIND/Debug/matlab');

% path to tools
tmpPATHS.tools = fullfile(tmpPATHS.code, 'tools/');
tmpPATHS.tools_nurbs = fullfile(tmpPATHS.code, 'tools/nurbs_toolbox/');

% path to DSC measurements
tmpPATHS.DSC204 = fullfile(tmpPATHS.root, 'DSC204_F1_Phoenix_Messungen/');
tmpPATHS.DSC204_measurements_pcm = fullfile(tmpPATHS.root, 'DSC204_F1_Phoenix_Messungen/Messungen/Messungen/');
tmpPATHS.DSC204_measurements_sap = fullfile(tmpPATHS.root, 'DSC204_F1_Phoenix_Messungen/Waermekapazitaet_Saphirmessung/');

% add everything to search path
for field = fieldnames(tmpPATHS)'
   fprintf('Adding path: %s', tmpPATHS.(field{1}));
   if exist(tmpPATHS.(field{1}), 'dir')
      addpath(tmpPATHS.(field{1}));
      fprintf(' ... OK!');
   else
      fprintf(' ... FAILED!');
   end
   fprintf('\n')
end
