function [spResPaths, contrResPaths] = extractMicrobeContributions_2026Gout(fbaFolder, saveDir)
% Efficiently extracts the microbe biomass shadow prices and relative
% abundances from a set of WBM FBA solutions and uses them to compute the
% potential microbial flux contributions, saving both to disk.
%
% USAGE:
%
%    [spResPaths, contrResPaths] = extractMicrobeContributions_2026Gout(fbaFolder, saveDir)
%
% INPUTS:
%    fbaFolder:      Char, directory containing the FBA solution `.mat` files
%                    (as written by testFluxExistingDiet_2026Gout / testFluxFoodItems_2026Gout).
%    saveDir:        Char, directory where the shadow-price and contribution
%                    results are written.
%
% OUTPUTS:
%    spResPaths:     Cell array of paths to the biomass shadow-price results,
%                    one `<rxn>.csv` per solved reaction.
%    contrResPaths:  Cell array of paths to the potential microbial
%                    contribution results, one `<rxn>.csv` per solved reaction.
%
% NOTE:
%
%    The reaction names are taken from the `rxns` field of the FBA solutions
%    (e.g. `173589_urate[bc]` or `DACHDiet_urate[bc]`), which matches the
%    column names of the fbaRes*.xlsx result tables. Germ-free WBM solutions
%    (file names containing `_gfWBM`) are excluded. Requires the COBRA Toolbox.
%
% .. Author: - Bram Nap

% Find paths to fba results
solPaths = fullfile(what(fbaFolder).path,what(fbaFolder).mat);
solPaths(contains(solPaths,'_gfWBM'))=[];

% Load fba solutions
solutions = cellfun( @(x) load(x,'ID','rxns','taxonNames','shadowPriceBIO','relAbundances'), solPaths);

% Find all microbial species
allSpecies = unique(vertcat(solutions(:).taxonNames));

% Reaction names, one per column of the solutions (diet/food x metabolite)
rxnNames = solutions(1).rxns(:);

% Preallocate result arrays
numRxns = numel(rxnNames);
numIDs = height(solutions);
numTaxa = length(allSpecies);

biomShadowPrices = nan(numIDs, numTaxa, numRxns);
relAbun = nan(numIDs, numTaxa);
microbeContributions = nan(numIDs, numTaxa, numRxns);

for i = 1:size(solutions, 1)

    % Find which species in the sample are present
    [index1, index2] = ismember(allSpecies, solutions(i).taxonNames);

    % Remove zero indices for non-present species in index2
    index2 = index2(index2 > 0);

    if ~isempty(index2)
        % Extract data
        samp_SP = solutions(i).shadowPriceBIO; % Extract shadow prices
        samp_ra = solutions(i).relAbundances; % Extract relative abundances

        % Filter microbes in model
        samp_SP = samp_SP(index2, :);
        samp_ra = samp_ra(index2);

        % Calculate microbe flux contribution potentials
        samp_c = samp_SP.*samp_ra;

        % Assign the shadow prices to the result array
        biomShadowPrices(i, index1, :) = samp_SP;
        relAbun(i,index1) = samp_ra;
        microbeContributions(i, index1, :) = samp_c * -1;
    end
end

% Obtain microbial biomass shadow prices for the predicted fluxes
createSpTab = @(x) array2table(biomShadowPrices(:,:,x),'RowNames', vertcat(solutions(:).ID), 'VariableNames',allSpecies);
shadowpriceTables = arrayfun(createSpTab,1:numRxns,'UniformOutput',false);

% Obtain total microbial component of the fluxes
createContrTab = @(x) array2table(microbeContributions(:,:,x),'RowNames', vertcat(solutions(:).ID), 'VariableNames',allSpecies);
contributionTables = arrayfun(createContrTab,1:numRxns,'UniformOutput',false);

% Create folders for shadow price and microbe contribution tables
spResFolder = fullfile(saveDir,'biomass_shadow_prices');  if ~isfolder(spResFolder); mkdir(spResFolder); end
contrResFolder = fullfile(saveDir,'potentialMicrobeContributions'); if ~isfolder(contrResFolder); mkdir(contrResFolder); end

% Save shadow price tables. Row names (model IDs) are written so that
% findMicrobialContributors_2026Gout can read them back with 'ReadRowNames'.
spResPaths = fullfile(spResFolder, append(rxnNames,'.csv'));
for i = 1:numRxns
    writetable(shadowpriceTables{i}, spResPaths{i}, 'WriteRowNames', true);
end

% Save microbe contribution tables
contrResPaths = fullfile(contrResFolder, append(rxnNames,'.csv'));
for i = 1:numRxns
    writetable(contributionTables{i}, contrResPaths{i}, 'WriteRowNames', true);
end
end
