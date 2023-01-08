% This file changes a rounded polytope in xml format into a polytope in mat
% format.
curFolder = fileparts(mfilename('fullpath'));
tmp = dir(fullfile(curFolder, '/PolyRound/output'));
matfiles = [];
for c = 1:size(tmp,1)
    if tmp(c).name(1) == 'p'
        matfiles = [matfiles tmp(c)];
    end
end
numPoly = 14;

for idx = 1:numPoly
    name = matfiles(idx).name;
    name = strcat('_', name, '_rounded.csv');
    model = struct;
    model.A = readmatrix(strcat(fullfile(curFolder, '/PolyRound/output/', matfiles(idx).name), '/A', name));
    model.b = readmatrix(strcat(fullfile(curFolder, '/PolyRound/output/', matfiles(idx).name), '/b', name));
    model.N = readmatrix(strcat(fullfile(curFolder, '/PolyRound/output/', matfiles(idx).name), '/N', name));
    model.p_shift = readmatrix(strcat(fullfile(curFolder, '/PolyRound/output/', matfiles(idx).name), '/p_shift', name));
    model.start = readmatrix(strcat(fullfile(curFolder, '/PolyRound/output/', matfiles(idx).name), '/start', name));

    parentFolder = fileparts(fileparts(curFolder));
    %save(fullfile(parentFolder, strcat('/Instances/1testbed-polyRound/', matfiles(idx).name, '.mat')), 'model')
end