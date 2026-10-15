function addUrateUptake2WBMs_2026Gout(modelDir, numWorkers)
% Adds the uptake of dietary urate from the gut into the blood to every
% whole-body metabolic model (WBM) in a folder and overwrites the model
% files. By default dietary urate only passes from the gut lumen to the
% feces; with these reactions it is absorbed in the small and large
% intestine, transported via the portal blood and liver into the
% circulation, and can be excreted in urine.
%
% USAGE:
%
%    addUrateUptake2WBMs_2026Gout(modelDir, numWorkers)
%
% INPUTS:
%    modelDir:      Char, folder with the WBMs (.mat), e.g. the foodDietWBMs
%                   folder written by setInSilicoDiet. The files are
%                   overwritten.
%
% OPTIONAL INPUTS:
%    numWorkers:    Numeric, number of parallel workers; 0 runs serially
%                   (default: 1).
%
% NOTE:
%
%    Uses addMetFromDiet2Urine2WBM_2026Gout(model, 'urate', 'Urate', param)
%    with param.Diet, param.SIabsorp and param.LIabsorp set to 'y'. Models
%    that already contain the small-intestine urate uptake reaction
%    (sIEC_URATEt[luI]) are skipped, so the function can be rerun safely.
%    Requires the COBRA Toolbox.
%
% .. Author: - Bram Nap

if nargin < 2 || isempty(numWorkers)
    numWorkers = 1;
end

wbmModels = what(modelDir).mat;
if isempty(wbmModels)
    error('addUrateUptake2WBMs_2026Gout:noModels', 'No .mat models found in %s', modelDir);
end

% Restore the COBRA environment on parallel workers (same pattern as
% testFluxFoodItems_2026Gout)
environment = getEnvironment();
if numWorkers > 0
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers);
    elseif poolobj.NumWorkers ~= numWorkers
        delete(poolobj);
        parpool(numWorkers);
    end
end

added = false(numel(wbmModels), 1);
parfor (i = 1:numel(wbmModels), numWorkers)
    restoreEnvironment(environment);
    fileName = fullfile(modelDir, wbmModels{i});
    model = load(fileName);

    % Skip models that already have the urate uptake route
    if any(strcmp(model.rxns, 'sIEC_URATEt[luI]'))
        continue
    end

    urateParam = struct('Diet', 'y', 'SIabsorp', 'y', 'LIabsorp', 'y');
    % addReactionsHH warns for each reaction that already exists (e.g.
    % Diet_EX_urate[d]); these are skipped and the warnings are not needed
    warnState = warning('off', 'all');
    model = addMetFromDiet2Urine2WBM_2026Gout(model, 'urate', 'Urate', urateParam);
    warning(warnState);

    saveModel(fileName, model);
    added(i) = true;
end
fprintf('Urate uptake added to %d of %d models in %s (%d already had it)\n', ...
    sum(added), numel(wbmModels), modelDir, sum(~added));
end

function saveModel(fileName, model)
% save() is not allowed directly inside parfor
save(fileName, '-struct', 'model');
end
