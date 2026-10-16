function [correlations, bestCorrPerSizeList, minbestCorrPerSizeList] = correlateClustersRestore_2026Gout(rxnsOfInterest, fluxData, microbeRelAbund, contributors)
% Function to create clusters from a set of microbes and correlate each
% cluster to the flux value for a given reaction
% How to use:
%   [correlations, bestCorrPerSizeList, minbestCorrPerSizeList] = correlateClustersRestore_2026Gout(rxnsOfInterest, fluxData, microbeRelAbund, contributors)
% Inputs:
%   rxnsOfInterest: Array, rxn IDs where clusters have to be found for
%   fluxData: Table, columns are rxn IDs, values are flux values
%   microbeRelAbund: Table, relative abundance used to create community
%                    microbiome models
%   contributors: Table, microbial contributors for each reaction as
%                 calculated from the shadow price in the
%                 findMicrobialContributorsRestore function
% Outputs:
%   correlations: Struct, each field is a rxn with the correlations results
%   bestCorrPerSizeList: Array, ordered by best Spearman's correlation
%                        coeffienct. For each cluster size, the best 
%                        combination of microbes is given. Positive
%                        correlation
%   minbestCorrPerSizeList: Array, ordered by best Spearman's correlation
%                        coeffienct. For each cluster size, the best 
%                        combination of microbes is given. Negative
%                        correlation
% Authors:
%   Tim Hensen: Nov 2025
%   Bram Nap:   Nov 2025

% Initialise structure to store all data in
correlations = struct();

% Initialise correlations arrays to store data in
bestCorrPerSizeList = cell(1,length(rxnsOfInterest));
minbestCorrPerSizeList = cell(1,length(rxnsOfInterest));

% For each reaction of interest
for i = 1:size(rxnsOfInterest,1)
    rxn = rxnsOfInterest{i};
    % Extract microbial species that influence that reaction
    if any(strcmpi(contributors.Properties.VariableNames, rxn))
        slice = contributors{contributors{:,rxn} ~= 0, 1};
        if any(strcmpi(slice, 'sex'))
            % Remove sex from the selection
            slice(strcmpi(slice,'sex')) = [];
        end

        % If there are less than 5 microbes, set maximum cluster size to the
        % number of microbes, otherwise set it to 5.
        if size(slice,1) < 5
            top = size(slice,1);
        else
            top = 5;
        end

        % Check if microbes influence the reaction
        if ~isempty(slice)
            % Create the clusters
            [clusterNames, clusterSize, clusterAbundances] = createClusters(top, slice, microbeRelAbund);
            clusterAbundances(isnan(clusterAbundances)) = 0;
            % Correlate the clusters
            clusterCorr = correlateClusters(fluxData, clusterAbundances, rxn, clusterNames, clusterSize);
            % Adjust the rxn name for valied fieldname
            rxnAdj = matlab.lang.makeValidName(rxn);
            % Add the results to the structure
            correlations.(rxnAdj) = clusterCorr;
            % Extract the the correlations and sort them into negative and postive
            % correlations
            [bestCorrPerSizeList{i},minbestCorrPerSizeList{i}] = extractTopClustersPerSize(clusterCorr,rxn);
        end
    end
end
end

function [clusterNames, clusterSize, clusterAbundances] = createClusters(top, topMicrobes, relAbun)
% Pre-calculate number of clusters for each size
clusterCounts = arrayfun(@(x) nchoosek(length(topMicrobes),x), 1:top);
numClusters = sum(clusterCounts);

% Pre-allocate arrays
clusterNames = strings(numClusters, top);
clusterAbundances = zeros(height(relAbun), numClusters, 'single');
clusterSize = zeros(numClusters, 1);

% Convert relAbun to numeric array once (outside loops)
numericAbun = single(table2array(relAbun(:, topMicrobes)));
numericAbun(isnan(numericAbun)) = 0;

% Initialize tracking variable
currentRow = 1;

% Generate combinations using a more efficient approach
for k = 1:top
    % Pre-allocate temporary arrays for batch processing
    combIdx = nchoosek(1:length(topMicrobes), k);
    batchSize = size(combIdx, 1);
    rowIndices = currentRow:(currentRow + batchSize - 1);

    % Set cluster sizes in batch
    clusterSize(rowIndices) = k;

    % Process names in batch of size k
    for j = 1:k
        clusterNames(rowIndices, j) = topMicrobes(combIdx(:, j));
    end

    % Create a sparse logical matrix for this batch
    batchIndicator = false(length(topMicrobes), batchSize);
    for i = 1:batchSize
        batchIndicator(combIdx(i,:), i) = true;
    end

    % Calculate abundances using matrix multiplication.
    % Note that for every matrix multiplication c(i,j) = a(i,k)*b(k,j) is defined as
    % for 1:k -> c(i,j) = sum( a(i,k) * b(k,j) ). This multiplication
    % speeds up the summation of the microbial abundances.
    clusterAbundances(:, rowIndices) = numericAbun * batchIndicator;

    % Update current row position
    currentRow = currentRow + batchSize;
end
end

function clusterCorr = correlateClusters(fluxValues, clusterAbundances, rxnName, clusterNames, clusterSize)
% Get reaction fluxes for metabolite
if any(contains(fluxValues.ID, 'gfWBM'))
    fluxValues(contains(fluxValues.ID, 'gfWBM'),:) = [];
end
rxnFlux = fluxValues.(rxnName);
disp(rxnName)
% Remove nan values in the fluxes and abundances for pairwise correlations
rowsToRemove = ~isnan(rxnFlux);
rxnFlux = rxnFlux(rowsToRemove);
clusterAbundances = clusterAbundances(rowsToRemove,:);

% Perform pearson correlation. Note that a spearman correlation would have
% been a bit more appropriate, however the spearman algorithm is much
% (much) slower and a linear regression of spearman and pearson results found that
% between 75 and 90% of the variance in the spearman results was captured
% in the pearson results. That is good enough for me.
RHO = corr(rxnFlux, clusterAbundances,'type','Spearman');

% Create table that combines the cluster names with the cluster sizes and
% the calculated correlations.
clusterCorr = array2table(clusterNames);
clusterCorr.Size = clusterSize;
clusterCorr.(rxnName) = RHO';

end

function [bestCorrPerSize,minbestCorrPerSize] = extractTopClustersPerSize(rxnClusterCorr,rxnName)
% Extract the top clusters for each cluster size

% Get the correlation coefficients
corrCoef = double(rxnClusterCorr.(rxnName));
% Find the group sizes
groupIds = rxnClusterCorr.Size;
% Get number of group sizes
numGroups = max(groupIds);

% Find the maximum and minimum correlation within each cluster size. Only
% the rows belonging to that size are searched (a sparse matrix would let
% its implicit zeros win when all correlations in a size share a sign).
maxIndices = nan(numGroups,1);
minIndices = nan(numGroups,1);
for g = 1:numGroups
    groupRows = find(groupIds == g);
    [maxVal, maxPos] = max(corrCoef(groupRows));
    [minVal, minPos] = min(corrCoef(groupRows));
    % Only keep clusters that correlate in the requested direction
    if maxVal > 0
        maxIndices(g) = groupRows(maxPos);
    end
    if minVal < 0
        minIndices(g) = groupRows(minPos);
    end
end
maxIndices = maxIndices(~isnan(maxIndices));
minIndices = minIndices(~isnan(minIndices));

% Extract the clusters with the highest correlations for each cluster size
bestCorrPerSize = rxnClusterCorr(maxIndices,:);
minbestCorrPerSize = rxnClusterCorr(minIndices,:);


% Perform the following steps to prepare the correlation tables
reVarsMicrobes = @(x) renamevars(x,x.Properties.VariableNames,replace(x.Properties.VariableNames,'clusterNames','Microbe_')); % Rename "cluster" to "Microbe"
reVar = @(x) renamevars(x,{'Size',rxnName},{'Cluster_size','Spearman_RHO'}); % Rename correlation column
addMet = @(x) addvars(x,repmat(cellstr(rxnName),height(x),1),'Before',1,'NewVariableNames','Metabolite'); % Specify associate metabolite as a separate column
addImpr = @(x) addvars(x,[nan(min(height(x),1),1); diff(abs(x.Spearman_RHO))],'After','Spearman_RHO','NewVariableNames','Absolute_RHO_improvement'); % Add the improved rho values for each increase in size
applyF = @(x) addImpr(addMet(reVar(reVarsMicrobes(x)))); % Combine defined table transformations

% Process top correlation tables
bestCorrPerSize = applyF(bestCorrPerSize);
minbestCorrPerSize = applyF(minbestCorrPerSize);

% With a single candidate microbe array2table names the column
% "clusterNames" without a number, so make it "Microbe_1"
if any(strcmp(bestCorrPerSize.Properties.VariableNames,'Microbe_'))
    bestCorrPerSize = renamevars(bestCorrPerSize,'Microbe_','Microbe_1');
    minbestCorrPerSize = renamevars(minbestCorrPerSize,'Microbe_','Microbe_1');
end

% Add correlation direction information
bestCorrPerSize = addvars(bestCorrPerSize,repmat({'Max'},height(bestCorrPerSize),1),'After','Metabolite','NewVariableNames','Direction');
minbestCorrPerSize = addvars(minbestCorrPerSize,repmat({'Min'},height(minbestCorrPerSize),1),'After','Metabolite','NewVariableNames','Direction');

% Add original cluster row index
bestCorrPerSize.clusterIndex = maxIndices;
minbestCorrPerSize.clusterIndex = minIndices;
end