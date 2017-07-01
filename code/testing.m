save_dir = '/home/argo/masterarbeit/simulationen-data/test/';

% check if directory where we want to save our data exist, if not create a
% folder
if exist(save_dir, 'dir') == 0
    disp('Creating folder!')
    mkdir(save_dir);
end

setup_path = strcat(save_dir, 'setup.mat');

fit_setup = struct()



return

save(setup_path, '-struct', 's1')



