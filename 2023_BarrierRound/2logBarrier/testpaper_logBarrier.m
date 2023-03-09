seed = 1;
rng(seed);

curFolder = fileparts(mfilename('fullpath'));
datapath = fullfile(fileparts(curFolder), '/Instances/2testbed-logBarrier', '*.mat');
matfiles = dir(datapath);
numPoly = size(matfiles, 1);
numSamples = 500;

list = 1;

for c = list
    % load polyoptes rounded by log-barrier
    model = load(fullfile(fileparts(curFolder), 'Instances/2testbed-logBarrier', matfiles(c).name));
    model = cellfun(@(x)(model.(x)), fieldnames(model));
    P = model;
    fprintf("Loading complete (%d/%d): %s\n", c, numPoly, matfiles(c).name(1:end-4))
    
    % sampling
    tic;
    samples = coordHAR(numSamples, size(P.A, 2)^2, P, [], 0.1);
    
    % return to original space
    samplesNew = P.x0 * ones(1, size(samples, 2));
    samplesNew(P.idx,:) = (samplesNew(P.idx,:) + P.p_shift) + P.N*samples;
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
    exps.roundTime = model.roundTime;
    save(fullfile(curFolder, strcat('/logBarrier_test/result_log_', matfiles(c).name)), 'exps')
end