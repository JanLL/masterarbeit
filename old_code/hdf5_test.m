filename = '/home/argo/Desktop/test.hdf5';

try
h5create(filename, '/Dset1', [10, 10]); 
catch
    fprintf('Dataset /Dset1 already exist!\n');
end

h5write(filename, '/Dset1', ones(10) .*3);

h5writeatt(filename, '/Dset1', 'N1', 1000);