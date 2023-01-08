seed = 1;
rng(seed);

curFolder = fileparts(mfilename('fullpath'));
datapath = fullfile(fileparts(curFolder), '/Instances/1testbed-polyRound', '*.mat');
matfiles = dir(datapath);
numPoly = size(matfiles, 1);
numSamples = 500;

list = 1:14;

for c = list
    % load polyoptes rounded by PolyRound
    model = load(fullfile(fileparts(curFolder), 'Instances/1testbed-polyRound', matfiles(c).name));
    model = cellfun(@(x)(model.(x)), fieldnames(model));
    P = model;
    fprintf("Loading complete (%d/%d): %s\n", c, numPoly, matfiles(c).name(1:end-4))
    
    % sampling
    tic;
    samples = coordHAR(numSamples, size(P.A, 2)^2, P, [], 0.1);
    
    % return to original space
    samplesNew = P.N*samples + P.p_shift;
    sampleTime = toc;
    
    % thinning to ensure independence
    ess = effectiveSampleSize(samplesNew);
    nSamples = min(ess);
    gap = round(size(samplesNew,2) / nSamples);
    samplesNew = samplesNew(:,gap:gap:end);
    
    % record results
    exps = struct; exps.dim = size(P.A, 2); 
    exps.ess = nSamples;
    exps.sampleTime = sampleTime;
    exps.step = size(samples, 2) * size(P.A, 2)^2;
    save(fullfile(curFolder, strcat('/PolyRound_test/result_PR_', matfiles(c).name)), 'exps')
end