seed = 1;
rng(seed);

curFolder = fileparts(mfilename('fullpath'));
datapath = fullfile(fileparts(curFolder), '/Instances/3testbed-isotropic', '*.mat');
matfiles = dir(datapath);
numPoly = size(matfiles, 1);
numSamples = 500;

list = 1:14;

for c = list
    % load polyoptes rounded by isotropic rounding
    model = load(fullfile(fileparts(curFolder), 'Instances/3testbed-isotropic', matfiles(c).name));
    model = cellfun(@(x)(model.(x)), fieldnames(model));
    P = model;
    fprintf("Loading complete (%d/%d): %s\n", c, numPoly, matfiles(c).name(1:end-4))
    
    % sampling
    tic;
    samples = coordHAR(numSamples, size(P.A, 2)^2, P, [], 0.1);

    % return to original space
    samplesNew = P.T*samples + P.y;
    sampleTime = toc;
    
    % check isotropy
    covmtx = checkCov(samples);
    e = eig(covmtx);
    %e = sort(e, 'descend');
    %plot(e, 'o')
    sprintf('max: %f, min: %f', max(e), min(e))

    % check feasibility
    %assert(min(samplesNew-model.lb, [], 'all') >= 0);
    %assert(min(model.ub-samplesNew, [], 'all') >= 0);
    %assert(min(vecnorm(model.S*samplesNew-model.b)) < 1e-8);
    
    % thinning to ensure independence
    ess = effectiveSampleSize(samplesNew);
    nSamples = min(ess);
    gap = round(size(samplesNew,2) / nSamples);
    samplesNew = samplesNew(:,gap:gap:end);
    
    % record results
    exps = struct; exps.dim = size(P.A, 2); 
    exps.ess = nSamples;
    exps.sampleTime = sampleTime;
    exps.eigs = e;
    exps.step = size(samples, 2) * size(P.A, 2)^2;
    exps.roundTime = model.roundTime;
    save(fullfile(curFolder, strcat('/isotropic_test/result_iso_', matfiles(c).name)), 'exps')
end

function covmtx = checkCov(samplesNew)
    ess = effectiveSampleSize(samplesNew);
    nSamples = min(ess);
    gap = round(size(samplesNew,2) / nSamples);
    samplesNew = samplesNew(:,gap:gap:end);

    o_shift = samplesNew - mean(samplesNew, 2); % now, mean is at the origin
    covmtx = o_shift*o_shift'/(nSamples-1);
end