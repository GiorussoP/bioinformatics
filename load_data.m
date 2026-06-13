% Load matrix A and vector b from a .mat file
if isfile('fibro_dataset.mat')
    load('fibro_dataset.mat', 'A', 'b', 'geneNames');
    fprintf('Loaded A, b and geneNames from fibro_dataset.mat successfully.\n');
else
    warning('File fibro_dataset.mat not found.');
end

% Load matrix Ar from a .mat file
if isfile('reduced_fibro_dataset.mat')
    load('reduced_fibro_dataset.mat', 'Ar');
    fprintf('Loaded Ar from reduced_fibro_dataset.mat successfully.\n');
else
    warning('File reduced_fibro_dataset.mat not found.');
end

