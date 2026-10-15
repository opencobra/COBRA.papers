function [maxCorr, minCorr] = findBestMicrobeClusters_2026Gout(fluxFilePath, microbeAbundancePath, enetResultPath, saveDir)
% Finds, for every predicted flux, the clusters (1 to 5 microbes) of the
% microbial contributors whose summed relative abundance correlates best
% (Spearman) with the flux, using correlateClustersRestore_2026Gout, and saves the
% best positive and negative cluster per cluster size with
% saveCorrTableRestore_2026Gout.
%
% USAGE:
%
%    [maxCorr, minCorr] = findBestMicrobeClusters_2026Gout(fluxFilePath, microbeAbundancePath, enetResultPath, saveDir)
%
% INPUTS:
%    fluxFilePath:          Char, flux table (fbaResDiet/fbaResFood xlsx); first
%                           column the model IDs, then one column per
%                           reaction, last column Sex.
%    microbeAbundancePath:  Char, relative abundance table of the microbes in
%                           each mWBM (first column the mWBM IDs).
%    enetResultPath:        Char, enetFluxMicrobeRes.csv written by
%                           findMicrobialContributors_2026Gout (long format: Reaction,
%                           Taxa, Beta, ...).
%    saveDir:               Char, directory where clusterCorrelations.xlsx is
%                           saved.
%
% OUTPUTS:
%    maxCorr:               Cell array, best positive correlation per cluster
%                           size for each reaction (see correlateClustersRestore_2026Gout).
%    minCorr:               Cell array, best negative correlation per cluster
%                           size for each reaction.
%
% NOTE:
%
%    correlateClustersRestore_2026Gout pairs fluxes and abundances by ROW POSITION, so
%    the abundance table is reduced to the models in the flux table, in the
%    same order. The model IDs are matched after removing quotes, 'foodItem'
%    and '.mat' from the flux-table IDs (as in findMicrobialContributors_2026Gout).
%    The elastic-net results are converted to the wide Beta table (taxa x
%    reactions, 0 when a taxon was not selected) that
%    correlateClustersRestore_2026Gout expects.
%
% .. Author: - Bram Nap

% Flux table with an ID column
fluxes = readtable(fluxFilePath, 'VariableNamingRule', 'preserve');
fluxes.Properties.VariableNames{1} = 'ID';
fluxes.ID = erase(fluxes.ID, {'''', '"', 'foodItem', '.mat'});

% Relative abundances in the same row order as the flux table
microbeRelAbund = readtable(microbeAbundancePath, 'VariableNamingRule', 'preserve');
microbeRelAbund.Properties.VariableNames{1} = 'ID';
[found, idx] = ismember(fluxes.ID, microbeRelAbund.ID);
if ~all(found)
    warning('findBestMicrobeClusters_2026Gout:missingAbundance', ...
        '%d model(s) not found in the abundance table and removed: %s', ...
        sum(~found), strjoin(fluxes.ID(~found)', ', '));
    fluxes = fluxes(found, :);
    idx = idx(found);
end
microbeRelAbund = microbeRelAbund(idx, :);

% Wide table of elastic-net coefficients: first column Taxa, one column per
% reaction, 0 when the taxon was not selected for the reaction
enet = readtable(enetResultPath, 'VariableNamingRule', 'preserve');
contributors = unstack(enet(:, {'Taxa', 'Reaction', 'Beta'}), 'Beta', 'Reaction', ...
    'VariableNamingRule', 'preserve');
contributors(:, 2:end) = fillmissing(contributors(:, 2:end), 'constant', 0);

% Reactions with microbial contributors
rxnsOfInterest = contributors.Properties.VariableNames(2:end)';

[~, maxCorr, minCorr] = correlateClustersRestore_2026Gout(rxnsOfInterest, fluxes, microbeRelAbund, contributors);

% Save (only when at least one cluster was found)
if ~isfolder(saveDir)
    mkdir(saveDir);
end
if any(cellfun(@(x) istable(x) && height(x) > 0, [maxCorr, minCorr]))
    saveCorrTableRestore_2026Gout(maxCorr, minCorr, saveDir);
else
    warning('findBestMicrobeClusters_2026Gout:noClusters', 'No microbial clusters found; nothing saved.');
end
end
