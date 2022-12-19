% Identify the most accurate model based on the results of driver_testiDNmodelPredictiveCapacity.m.
% The information gathered from each model is used to identify two sets of conditions:
%   - The conditions for the model with the highest score;
%   - The conditions with the highest score on each condition used.
%
% The most accurate models are saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1

clear

%param.approach ='SD';
%param.approach ='mean';
param.approach ='UptSec';

legends = 0;

allPlots = 0;

letters = 'A':'Z';

format long g

writeSBML = 0;

% Define results directory
currentDir = pwd;
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

%define the directory with the ensemble of models
modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'multidimensionalModelGeneration'];

% Load multidimensionalComparisonStats
load([pathSave filesep 'multidimensionalComparisonStats' param.approach '.mat'])
directoriesWithModels = unique(multidimensionalComparisonStats.dirName);

[~, ia] = sort(multidimensionalComparisonStats.quantitativeBoth, 'ascend');
multidimensionalComparisonStats = addvars(multidimensionalComparisonStats,zeros(size(multidimensionalComparisonStats,1),1),'NewVariableNames','quantitativeBothScore','After','quantitativeBoth');
for i=1:size(multidimensionalComparisonStats,1)
    multidimensionalComparisonStats.quantitativeBothScore(ia(i))=1 - i/size(multidimensionalComparisonStats,1);
end

[~, ia] = sort(multidimensionalComparisonStats.quantitativeBoth, 'descend');
multidimensionalComparisonStats = addvars(multidimensionalComparisonStats,zeros(size(multidimensionalComparisonStats,1),1),'NewVariableNames','quantitativeBothScore2','After','quantitativeBoth');
for i=1:size(multidimensionalComparisonStats,1)
    multidimensionalComparisonStats.quantitativeBothScore2(ia(i))= i/size(multidimensionalComparisonStats,1);
end

multidimensionalComparisonStats.rowRankDeficiency = zeros(size(multidimensionalComparisonStats,1),1);
multidimensionalComparisonStats = addvars(multidimensionalComparisonStats,zeros(size(multidimensionalComparisonStats,1),1),'NewVariableNames','selectionScore','After','objective');

%% Score for selection of best model
for i=1:size(multidimensionalComparisonStats,1)
    multidimensionalComparisonStats.rowRankDeficiency(i)=multidimensionalComparisonStats.nOfmets(i) - multidimensionalComparisonStats.rankOfS(i);
    multidimensionalComparisonStats.selectionScore(i) = sum(multidimensionalComparisonStats.qualitativeBoth(i)  + multidimensionalComparisonStats.spearmanBoth(i) + multidimensionalComparisonStats.quantitativeBothScore(i))/3;
    multidimensionalComparisonStats.selectionScore2(i) = sum(multidimensionalComparisonStats.qualitativeBoth(i)  + multidimensionalComparisonStats.spearmanBoth(i) + multidimensionalComparisonStats.quantitativeBothScore2(i))/3;
end
multidimensionalComparisonStats = sortrows(multidimensionalComparisonStats,'selectionScore','descend');

%% 3D plot of scores
if 0
    figure
    plot3(multidimensionalComparisonStats.qualitativeBoth,log10(multidimensionalComparisonStats.quantitativeBoth),multidimensionalComparisonStats.spearmanBoth,'.')
    xlabel('Correct/total')
    ylabel('Euclidean distance')
    zlabel('Spearman rho')
    return
end

%% Determine the models with the greatest predictive capacity across all objectives
models = unique(multidimensionalComparisonStats.dirName);

[qualitativeAccuracy, quantitativeAccuracy, spearmanAccuracy,selectionScore] = deal(zeros(size(models)));
for i = 1:length(models)
    
    % Find rows
    rowsInComparison = strcmp(multidimensionalComparisonStats.dirName, ...
        models{i});
    
    % Qualitative accuracy
    qualitativeScoreBoth = multidimensionalComparisonStats.qualitativeBoth(rowsInComparison);
    qualitativeAccuracy(i) = mean(qualitativeScoreBoth,'omitnan');
    
    % Euclidean distance score
    quantitativeScoreBoth = multidimensionalComparisonStats.quantitativeBothScore(rowsInComparison);
    quantitativeAccuracy(i) = mean(quantitativeScoreBoth,'omitnan');
    
    % Spearman rho
    spearmanScoreBoth = multidimensionalComparisonStats.spearmanBoth(rowsInComparison);
    spearmanAccuracy(i) = mean(spearmanScoreBoth,'omitnan');
    
    % Selection score
    selectionScoreBoth = multidimensionalComparisonStats.selectionScore(rowsInComparison);
    % Ignores outlier, and NaN values
    toIgnore = isnan(selectionScoreBoth);
    selectionScore(i,1) = mean(selectionScoreBoth,'omitnan');
    
end

% best nlt models
nlt =10;
qualitativeAccuracy(isnan(qualitativeAccuracy)) = 0;
[~, ia] = sort(qualitativeAccuracy, 'descend');
disp(table(models(ia(1:nlt)),qualitativeAccuracy(ia(1:nlt)),'VariableNames',{'model','Quantitative accuracy'}))

spearmanAccuracy(isnan(spearmanAccuracy)) = 0;
[~, ia] = sort(spearmanAccuracy, 'descend');
disp(table(models(ia(1:nlt)),spearmanAccuracy(ia(1:nlt)),'VariableNames',{'model','Spearman rho'}))

quantitativeAccuracy(isnan(quantitativeAccuracy)) = 0;
[~, ia] = sort(quantitativeAccuracy, 'ascend');
disp(table(models(ia(1:nlt)),quantitativeAccuracy(ia(1:nlt)),'VariableNames',{'model','Euclidean distance'}))

selectionScore(isnan(selectionScore)) = 0;
[~, ia] = sort(selectionScore, 'descend');
disp(table(models(ia(1:nlt)),selectionScore(ia(1:nlt)),'VariableNames',{'model','Selection score'}))

%% Determine the objectives with the greatest predictive fidelity across all models
objectives = unique(multidimensionalComparisonStats.objective);

[qualitativeAccuracy, quantitativeAccuracy, spearmanAccuracy, selectionScore] = deal(zeros(size(objectives)));
for i = 1:length(objectives)
    
    % Find rows
    rowsInComparison = strcmp(multidimensionalComparisonStats.objective, ...
    objectives{i});
    
    % Qualitative accuracy
    qualitativeScoreBoth = multidimensionalComparisonStats.qualitativeBoth(rowsInComparison);
    qualitativeAccuracy(i) = mean(qualitativeScoreBoth,'omitnan');
    
    % Euclidean distance
    quantitativeScoreBoth = multidimensionalComparisonStats.quantitativeBothScore(rowsInComparison);
    quantitativeAccuracy(i) = mean(quantitativeScoreBoth,'omitnan');
    
    % Spearman rho
    spearmanScoreBoth = multidimensionalComparisonStats.spearmanBoth(rowsInComparison);
    spearmanAccuracy(i) = mean(spearmanScoreBoth,'omitnan');
    
    % Selection score
    selectionScoreBoth = multidimensionalComparisonStats.selectionScore(rowsInComparison);
    selectionScore(i) = mean(selectionScoreBoth,'omitnan');
end

% best objectives
qualitativeAccuracy(isnan(qualitativeAccuracy)) = 0;
[~, ia] = sort(qualitativeAccuracy, 'descend');
disp(table(objectives(ia),qualitativeAccuracy(ia),'VariableNames',{'objective','Qualitative_accuracy'}))

spearmanAccuracy(isnan(spearmanAccuracy)) = 0;
[~, ia] = sort(spearmanAccuracy, 'descend');
disp(table(objectives(ia),spearmanAccuracy(ia),'VariableNames',{'objective','Spearman_rho'}))

quantitativeAccuracy(isnan(quantitativeAccuracy)) = 0;
[~, ia] = sort(quantitativeAccuracy, 'ascend');
disp(table(objectives(ia),quantitativeAccuracy(ia),'VariableNames',{'objective','Euclidean_score'}))

selectionScore(isnan(selectionScore)) = 0;
[~, ia] = sort(selectionScore, 'descend');
disp(table(objectives(ia),selectionScore(ia),'VariableNames',{'objective','Selection_score'}))

% make objectives table
T = table(objectives, qualitativeAccuracy, spearmanAccuracy, quantitativeAccuracy,selectionScore);
T = sortrows(T,'selectionScore','descend');
table2latex(T, [pathSave filesep 'objectives.tex']);

%save([pathSave filesep 'topObjectives' param.approach], 'topObjectives')

%% Subset of models+objectives to plot the effects of varying extraction parameters
if 0
    % Models with more than 1000 reactions
    modelsOverMinNumRxns = multidimensionalComparisonStats.rankOfS > 1000;
else
    modelsOverMinNumRxns = multidimensionalComparisonStats.rankOfS >= min(multidimensionalComparisonStats.rankOfS);
end

%% Select one or a set of objectives to evaluate
%objectivesBool = ismember(multidimensionalComparisonStats.objective, topObjectives(1));
%objectivesBool = ismember(multidimensionalComparisonStats.objective, topObjectives(1:3));
%objectivesBool = ismember(multidimensionalComparisonStats.objective, topObjectives);
%objectivesBool = ismember(multidimensionalComparisonStats.objective, objectives(ismember(objectives,topObjectives)));
%objectivesBool = ismember(multidimensionalComparisonStats.objective, objectives); %all objectives
objectivesBool = ismember(multidimensionalComparisonStats.objective, 'unWeightedTCBMflux');

[totalQualitativeScoreModels, totalSpearmanScoreModels, totalEuclideanScoreModels, noOfRxns, noOfMets, ...
    rankOfS, rowRankDeficiency, spearmanScoreModels,qualitativeScoreModels,euclideanScoreModels] = deal(zeros(size(directoriesWithModels)));
for i = 1:length(directoriesWithModels)
    
    dirBool = strcmp(multidimensionalComparisonStats.dirName, directoriesWithModels{i});
    rowsInComparison = dirBool & objectivesBool & modelsOverMinNumRxns;
    
    if any(rowsInComparison)
        
        % Qualitative score
        qualitativeScoreBoth = multidimensionalComparisonStats.qualitativeBoth(rowsInComparison);
        qualitativeScoreModelSec = multidimensionalComparisonStats.qualitativeModelSec(rowsInComparison);
        qualitativeScoreModelUpt = multidimensionalComparisonStats.qualitativeModelUpt(rowsInComparison);
        % Ignores outlier, and NaN values
        %toIgnore = isoutlier(qualitativeScoreBoth) | isnan(qualitativeScoreModelSec) | isnan(qualitativeScoreModelUpt);
        toIgnore = isnan(qualitativeScoreModelSec) | isnan(qualitativeScoreModelUpt);
        qualitativeScoreModels(i, 1) = mean(qualitativeScoreBoth(~toIgnore));
        
        % Quantitative score
        spearmanScoreBoth = multidimensionalComparisonStats.spearmanBoth(rowsInComparison);
        spearmanScoreModelSec = multidimensionalComparisonStats.spearmanModelSec(rowsInComparison);
        spearmanScoreModelUpt = multidimensionalComparisonStats.spearmanModelUpt(rowsInComparison);
        % Ignores outlier, and NaN values
        %toIgnore = isoutlier(spearmanScoreBoth) | isnan(spearmanScoreModelSec) | isnan(spearmanScoreModelUpt);
        toIgnore = isnan(spearmanScoreModelSec) | isnan(spearmanScoreModelUpt);
        spearmanScoreModels(i, 1) = mean(spearmanScoreBoth(~toIgnore));
        
        % Quantitative score
        euclideanScoreBoth = multidimensionalComparisonStats.quantitativeBoth(rowsInComparison);
        euclideanScoreModelSec = multidimensionalComparisonStats.quantitativeModelSec(rowsInComparison);
        euclideanScoreModelUpt = multidimensionalComparisonStats.quantitativeModelUpt(rowsInComparison);
        % Ignores outlier, and NaN values
        %toIgnore = isoutlier(euclideanScoreBoth) | isnan(euclideanScoreModelSec) | isnan(euclideanScoreModelUpt);
        toIgnore = isnan(euclideanScoreModelSec) | isnan(euclideanScoreModelUpt);
        euclideanScoreModels(i, 1) = mean(euclideanScoreBoth(~toIgnore));
        
        % Number of reactions and metabolites
        noOfMets(i) = unique(multidimensionalComparisonStats.nOfmets(dirBool));
        noOfRxns(i) = unique(multidimensionalComparisonStats.nOfrxns(dirBool));
        rankOfS(i) = unique(multidimensionalComparisonStats.rankOfS(dirBool));
        rowRankDeficiency(i) = unique(multidimensionalComparisonStats.rowRankDeficiency(dirBool));
    else
        
        qualitativeScoreModels(i, 1) = 0;
        spearmanScoreModels(i, 1) = 0;
        euclideanScoreModels(i, 1) = 0;
        % Number of reactions and metabolites
        noOfMets(i) = 0;
        noOfRxns(i) = 0;
        rankOfS(i) = 0;
        rowRankDeficiency(i) = 0;
        
    end
    totalQualitativeScoreModels(i,1) = mean(nonzeros(qualitativeScoreModels(i, :)));
    totalSpearmanScoreModels(i,1) = mean(nonzeros(spearmanScoreModels(i, :)));
    totalEuclideanScoreModels(i,1) = mean(nonzeros(euclideanScoreModels(i, :)));
end

% Set zeros to NaN
qualitativeScoreModels(qualitativeScoreModels == 0) = NaN;
spearmanScoreModels(spearmanScoreModels == 0) = NaN;
noOfMets(noOfMets == 0) = NaN;
noOfRxns(noOfRxns == 0) = NaN;
rankOfS(rankOfS == 0) = NaN;
rowRankDeficiency(rowRankDeficiency == 0) = NaN;


%% mean difference for each condition
% Separate each condition in the model
directoriesComparison = cell(size(directoriesWithModels, 1), ...
    length(split(directoriesWithModels{1}, '_')));
for i = 1:size(directoriesWithModels, 1)
    directoriesComparison(i, :) = split(directoriesWithModels{i}, '_')';
end

% Identification of data names and groups
c = 0;
fields = cell(size(unique(directoriesComparison), 1), 1);
for i = 1:size(directoriesComparison, 2)
    currentConditions = unique(directoriesComparison(:, i));
    for j = 1:size(currentConditions, 1)
        c = c + 1;
        fields{c, 1} = [currentConditions{j} '_' num2str(i)];
    end
end
dataGroups = regexp(fields, '(\d+)(?!.*\d)', 'match');
dataGroups = [dataGroups{:}]';

% Adjustments to the labels
expression = 'limit|GenesT|Ions|media|Relaxed|transcriptomics';
plotLabels = regexprep(fields, expression, '');
plotLabels = regexprep(plotLabels, 'deleteModelGenes', 'all-rxn');
plotLabels = regexprep(plotLabels, 'oneRxnPerActiveGene', '1-rxn');
plotLabels = regexprep(plotLabels, 'T', '');
plotLabels = regexprep(plotLabels, 'Boundary.100000', '1e5');
plotLabels = regexprep(plotLabels, 'Boundary.10000', '1e4');
plotLabels = regexprep(plotLabels, 'Boundary.1000', '1e3');
plotLabels = regexprep(plotLabels, 'inactive', 'yes');
plotLabels = regexprep(plotLabels, 'NoInactive', 'no');
plotLabels = regexprep(plotLabels, 'Over([^;]*).*?(?=\_)', '');


figure
if allPlots
    t = tiledlayout(3,4);
else
    t = tiledlayout(2,3);
end
t.TileSpacing = 'compact';
k=1;
%% Qualitative accuracy
nexttile
hold on
limitAx = 0;
for i = 1:length(unique(dataGroups))
    groups{i} = ['group' num2str(i)];
end

% Plot the labels according the mean difference
conditionsAcc = zeros(size(plotLabels));
for i = 1:size(plotLabels, 1)
    parameterLabel = regexprep(fields{i}, '\_(.*)', '');
    groupBool = sum(ismember(directoriesComparison, parameterLabel), 2);
    nanBool = isnan(totalQualitativeScoreModels);
    meanQualitative = (mean(totalQualitativeScoreModels(groupBool & ~nanBool)) - ...
        mean(totalQualitativeScoreModels(~nanBool)));
    h1 = text(meanQualitative, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
        '\_(.*)', '')], 'FontSize', 12,'Color', [0 0.4510 0.7412]);
    set(h1, 'Rotation', 45)
    if limitAx < abs(meanQualitative)
        limitAx = abs(meanQualitative);
    end
    conditionsAcc(i) = meanQualitative;
    if i == size(plotLabels, 1)
        text(0, 0.5, ['  Mean: ' num2str(round(mean(totalQualitativeScoreModels(~isnan(...
            totalQualitativeScoreModels))), 2))], 'FontSize', 14);
    end
end
xline(0, '--r');
axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
hold off
yticks(1:7)
conditionLabels = {'Extraction algorithm', ...
    '#Reactions per active gene', 'Transcript. threshold, log_2(FPKM)', ...
    'Max. default flux (uMol/gDW/hr)', 'Transcript. inactive genes', 'Ion exchange',...
    'Data preference'};
yticklabels(conditionLabels)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14, 'fontweight', 'bold')
l=1;
title([letters(l) '. Correct/total predictions'], 'FontSize', 14)
l=l+1;
set(gca,'fontsize',14,'FontWeight','bold')
xlabel('Mean difference', 'FontSize', 12, 'fontweight', 'bold')
grid minor
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';


if allPlots || 0
    %% Euclidean distance subplot
    nexttile
    hold on
    limitAx = 0;
    conditionsED = zeros(size(plotLabels));
    for i = 1:size(plotLabels, 1)
        % Plot the labels according the mean difference
        parameterLabel = regexprep(fields{i}, '\_(.*)', '');
        groupBool = sum(ismember(directoriesComparison, parameterLabel), 2);
        nanBool = isnan(totalEuclideanScoreModels);
        meanQuantitative = (mean(totalEuclideanScoreModels(groupBool & ~nanBool)) - ...
            mean(totalEuclideanScoreModels(~nanBool)));
        h1 = text(meanQuantitative, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
            '\_(.*)', '')], 'FontSize', 12,'Color', [0.8510 0.3294 0.1020]);
        set(h1, 'Rotation', 45)
        if limitAx < abs(meanQuantitative)
            limitAx = abs(meanQuantitative);
        end
        conditionsED(i) = meanQuantitative;
        if i == size(plotLabels, 1)
            text(0, 0.5, ['  Mean: ' num2str(round(mean(totalEuclideanScoreModels(~nanBool)), 2))...
                ], 'FontSize', 14);
        end
    end
    xline(0, '--r');
    axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
    hold off
    title([letters(l) '. Weighted Euclidean distance'], 'FontSize', 14)
    l=l+1;
    xlabel('Mean difference', 'FontSize', 12, 'fontweight', 'bold')
    yticks(1:7)
    yticklabels({'','', '', '', '', '', ''})
    grid minor
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
end

if allPlots || 1
    %% Spearman rho subplot
    nexttile
    hold on
    limitAx = 0;
    conditionsED = zeros(size(plotLabels));
    for i = 1:size(plotLabels, 1)
        % Plot the labels according the mean difference
        parameterLabel = regexprep(fields{i}, '\_(.*)', '');
        groupBool = sum(ismember(directoriesComparison, parameterLabel), 2);
        nanBool = isnan(totalSpearmanScoreModels);
        meanQuantitative = (mean(totalSpearmanScoreModels(groupBool & ~nanBool)) - ...
            mean(totalSpearmanScoreModels(~nanBool)));
        h1 = text(meanQuantitative, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
            '\_(.*)', '')], 'FontSize', 12,'Color', [0.8510 0.3294 0.1020]);
        set(h1, 'Rotation', 45)
        if limitAx < abs(meanQuantitative)
            limitAx = abs(meanQuantitative);
        end
        conditionsED(i) = meanQuantitative;
        if i == size(plotLabels, 1)
            text(0, 0.5, ['  Mean: ' num2str(round(mean(totalSpearmanScoreModels(~nanBool)), 2))...
                ], 'FontSize', 14);
        end
    end
    xline(0, '--r');
    axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
    hold off
    title([letters(l) '. Spearman rho'], 'FontSize', 14)
    l=l+1;
    set(gca,'fontsize',14,'FontWeight','bold')
    xlabel('Mean difference', 'FontSize', 14, 'fontweight', 'bold')
    yticks(1:7)
    yticklabels({'','', '', '', '', '', ''})
    grid minor
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
end

if allPlots || 1
    %% rank subplot
    nexttile
    hold on
    limitAx = 0;
    conditionsRank = zeros(size(plotLabels));
    for i = 1:size(plotLabels, 1)
        % Plot the labels according the mean difference
        parameterLabel = regexprep(fields{i}, '\_(.*)', '');
        [idxConditionSpecific, ~] = find(ismember(directoriesComparison, parameterLabel));
        nanBool = isnan(totalSpearmanScoreModels);
        meanQuantitative = nanmean(rankOfS(idxConditionSpecific)) - nanmean(rankOfS);
        h1 = text(meanQuantitative, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
            '\_(.*)', '')], 'FontSize', 14,'Color', [0.9290, 0.6940, 0.1250]);
        set(h1, 'Rotation', 45)
        if limitAx < abs(meanQuantitative)
            limitAx = abs(meanQuantitative);
        end
        conditionsRank(i) = meanQuantitative;
        if i == size(plotLabels, 1)
            text(0, 0.5, ['  Mean: ' num2str(round(nanmean(rankOfS)))...
                ], 'FontSize', 14);
        end
    end
    xline(0, '--r');
    axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
    hold off
    title([letters(l) '. Stoichiometric matrix rank'], 'FontSize', 14)
    l=l+1;
    set(gca,'fontsize',14,'FontWeight','bold')
    xlabel('Mean difference', 'FontSize', 14, 'fontweight', 'bold')
    yticks(1:7)
    yticklabels({'', '', '', '', '', '', ''})
    grid minor
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
end

if allPlots
    %% row rank deficiency subplot
    nexttile
    hold on
    limitAx = 0;
    conditionsRank = zeros(size(plotLabels));
    for i = 1:size(plotLabels, 1)
        % Plot the labels according the mean difference
        parameterLabel = regexprep(fields{i}, '\_(.*)', '');
        [idxConditionSpecific, ~] = find(ismember(directoriesComparison, parameterLabel));
        nanBool = isnan(totalSpearmanScoreModels);
        meanQuantitative = nanmean(rowRankDeficiency(idxConditionSpecific)) - nanmean(rowRankDeficiency);
        h1 = text(meanQuantitative, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
            '\_(.*)', '')], 'FontSize', 14,'Color', [0.9290, 0.6940, 0.1250]);
        set(h1, 'Rotation', 45)
        if limitAx < abs(meanQuantitative)
            limitAx = abs(meanQuantitative);
        end
        conditionsRank(i) = meanQuantitative;
        if i == size(plotLabels, 1)
            text(0, 0.5, ['  Mean: ' num2str(round(nanmean(rowRankDeficiency)))...
                ], 'FontSize', 14);
        end
    end
    xline(0, '--r');
    axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
    hold off
    title('D. Row rank deficiency', 'FontSize', 14)
    xlabel('Mean difference', 'FontSize', 14, 'fontweight', 'bold')
    yticks(1:7)
    yticklabels({'', '', '', '', '', '', ''})
    grid minor
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
end

% savefig([pathSave filesep 'conditions'])
% saveas(gcf,[pathSave filesep 'conditions'],'png')


%% iDopaNeuroCT

% Identify cell-type-specific models
cellTypeSpecificBool = contains(multidimensionalComparisonStats.preferenceCurationOrOmics, ...
    'curationOverOmics');
% Identify the max predictive capacity in cell-type-specific models
maxCTpredictiveCapacity = max(multidimensionalComparisonStats.selectionScore(...
    cellTypeSpecificBool & modelsOverMinNumRxns & objectivesBool));
% Identify the directory of the most accurate models
iDopaNeuroCTBool = multidimensionalComparisonStats.selectionScore == maxCTpredictiveCapacity & ...
    cellTypeSpecificBool & modelsOverMinNumRxns & objectivesBool;
if nnz(iDopaNeuroCTBool)>1
    iDopaNeuroCTBool = multidimensionalComparisonStats.selectionScore == maxCTpredictiveCapacity & ...
        cellTypeSpecificBool & modelsOverMinNumRxns & objectivesBool & multidimensionalComparisonStats.ATPtm;
end
if nnz(iDopaNeuroCTBool)==0
    error('could not identify iDopaNeuroCTBool')
end


iDopaNeuroCTDir = multidimensionalComparisonStats.dirName(iDopaNeuroCTBool);
% Predictive capacity
iDopaNeuroCTData(1, 1) = unique(multidimensionalComparisonStats.nOfmets(iDopaNeuroCTBool));
iDopaNeuroCTData(2, 1) = unique(multidimensionalComparisonStats.nOfrxns(iDopaNeuroCTBool));
iDopaNeuroCTData(3, 1) = unique(multidimensionalComparisonStats.rankOfS(iDopaNeuroCTBool));
iDopaNeuroCTData(4, 1) = multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCTBool);
iDopaNeuroCTData(5, 1) = multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCTBool);
iDopaNeuroCTData(6, 1) = multidimensionalComparisonStats.quantitativeBoth(iDopaNeuroCTBool);

% Coordinates figures E and F
iDNCT_Ecoordinates = [iDopaNeuroCTData(4, 1) iDopaNeuroCTData(5, 1)];
iDNCT_Fcoordinates = [iDopaNeuroCTData(4, 1) iDopaNeuroCTData(3, 1)];

% Table
iDN1ModelOFValidation = [round(multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCTBool), 2)'; ...
    round(multidimensionalComparisonStats.qualitativeModelSec(iDopaNeuroCTBool), 2)'; ...
    round(multidimensionalComparisonStats.qualitativeModelUpt(iDopaNeuroCTBool), 2)'; ...
    round(multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCTBool), 2)'; ...
    round(multidimensionalComparisonStats.spearmanModelSec(iDopaNeuroCTBool), 2)'; ...
    round(multidimensionalComparisonStats.spearmanModelUpt(iDopaNeuroCTBool), 2)'];
iDopaNeuroCTValidationTable = table(iDN1ModelOFValidation(:, 1), ...
    'VariableNames', ...
    multidimensionalComparisonStats.objective(iDopaNeuroCTBool),...
    'RowNames',...
    {'qualitativeBoth'; ...
    'qualitativeModelSec';...
    'qualitativeModelUpt';...
    'spearmanBoth';...
    'spearmanModelSec';...
    'spearmanModelUpt'});
disp(iDopaNeuroCTValidationTable)

% Copy to results folder

if ~isfolder([pathSave filesep 'iDopaNeuroCT'])
    mkdir([pathSave filesep 'iDopaNeuroCT'])
    copyfile([modelsDir filesep char(iDopaNeuroCTDir)], [pathSave filesep 'iDopaNeuroCT'])
    load([pathSave filesep 'iDopaNeuroCT' filesep 'Model'])
    iDopaNeuroCT = Model;
    save([pathSave filesep 'iDopaNeuroCT' filesep 'iDopaNeuroCT'], 'iDopaNeuroCT')
    delete([pathSave filesep 'iDopaNeuroCT' filesep 'Model.mat'])
    
    if writeSBML
        cd([pathSave filesep 'iDopaNeuroCT'])
        writeCbModel(iDopaNeuroCT, 'format','sbml', 'fileName', 'iDopaNeuroCT')
        writeCbModel(iDopaNeuroCT,'xls','iDopaNeuroCT')
        cd(currentDir)
        
    end
end
load([pathSave filesep 'iDopaNeuroCT' filesep 'iDopaNeuroCT' '.mat'])
eval(['model = ' 'iDopaNeuroCT' ';'])
iDopaNeuroCTStats = modelStatistics(model);
iDopaNeuroCTData(7, 1) = iDopaNeuroCTStats.Statistics(strcmp(iDopaNeuroCTStats.Properties.RowNames,'Internal reactions'));
iDopaNeuroCTData(8, 1) = iDopaNeuroCTStats.Statistics(strcmp(iDopaNeuroCTStats.Properties.RowNames,'External reactions'));
iDopaNeuroCTData(9, 1) = iDopaNeuroCTStats.Statistics(strcmp(iDopaNeuroCTStats.Properties.RowNames,'Genes'));

% Print conditions
fid2 = fopen([pathSave filesep 'iDopaNeuroCT' filesep 'conditions.txt'], 'w');
fprintf(fid2, '%s\n', iDopaNeuroCTDir);
fclose(fid2);

%% iDopaNeuroC

% Identify condition-specific models
conditionSpecificBool = contains(multidimensionalComparisonStats.preferenceCurationOrOmics, 'omicsOverCuration');
% Identify the max predictive capacity in condition-specific models
maxCpredictiveCapacity = max(multidimensionalComparisonStats.selectionScore(...
    conditionSpecificBool & modelsOverMinNumRxns & objectivesBool));
% Identify the directory of the most accurate models
iDopaNeuroCIdx = find(...
    multidimensionalComparisonStats.selectionScore == maxCpredictiveCapacity & ...
    conditionSpecificBool & modelsOverMinNumRxns & objectivesBool);
if length(iDopaNeuroCIdx) > 1
    maxQuantitative = find(multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCIdx) == ...
        max(multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCIdx)));
    iDopaNeuroCIdx = iDopaNeuroCIdx(maxQuantitative);
end
iDopaNeuroCDir = multidimensionalComparisonStats.dirName(iDopaNeuroCIdx);
% Predictive capacity
iDopaNeuroCData(1, 1) = unique(multidimensionalComparisonStats.nOfmets(iDopaNeuroCIdx));
iDopaNeuroCData(2, 1) = unique(multidimensionalComparisonStats.nOfrxns(iDopaNeuroCIdx));
iDopaNeuroCData(3, 1) = unique(multidimensionalComparisonStats.rankOfS(iDopaNeuroCIdx));
iDopaNeuroCData(4, 1) = multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCIdx);
iDopaNeuroCData(5, 1) = multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCIdx);
iDopaNeuroCData(6, 1) = multidimensionalComparisonStats.quantitativeBoth(iDopaNeuroCIdx);
% Coordinates figures E and F
iDNC_Ecoordinates = [iDopaNeuroCData(4, 1) iDopaNeuroCData(5, 1)];
iDNC_Fcoordinates = [iDopaNeuroCData(4, 1) iDopaNeuroCData(3, 1)];

% Table
iDN1ModelOFValidation = [round(multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCIdx), 2)'; ...
    round(multidimensionalComparisonStats.qualitativeModelSec(iDopaNeuroCIdx), 2)'; ...
    round(multidimensionalComparisonStats.qualitativeModelUpt(iDopaNeuroCIdx), 2)'; ...
    round(multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCIdx), 2)'; ...
    round(multidimensionalComparisonStats.spearmanModelSec(iDopaNeuroCIdx), 2)'; ...
    round(multidimensionalComparisonStats.spearmanModelUpt(iDopaNeuroCIdx), 2)'];
iDopaNeuroCValidationTable = table(iDN1ModelOFValidation(:, 1), ...
    'VariableNames', ...
    multidimensionalComparisonStats.objective(iDopaNeuroCIdx),...
    'RowNames',...
    {'qualitativeBoth'; ...
    'qualitativeModelSec';...
    'qualitativeModelUpt';...
    'spearmanBoth';...
    'spearmanModelSec';...
    'spearmanModelUpt'});
disp(iDopaNeuroCValidationTable)

% Copy to results folder
if ~isfolder([pathSave filesep 'iDopaNeuroC'])
    mkdir([pathSave filesep 'iDopaNeuroC'])
    
    copyfile([modelsDir filesep char(iDopaNeuroCDir)], [pathSave filesep 'iDopaNeuroC'])
    load([pathSave filesep 'iDopaNeuroC' filesep 'Model'])
    iDopaNeuroC = Model;
    save([pathSave filesep 'iDopaNeuroC' filesep 'iDopaNeuroC'], 'iDopaNeuroC')
    delete([pathSave filesep 'iDopaNeuroC' filesep 'Model.mat'])
    
    if writeSBML
        cd([pathSave filesep 'iDopaNeuroC'])
        writeCbModel(iDopaNeuroC, 'format','sbml', 'fileName', 'iDopaNeuroC')
        writeCbModel(iDopaNeuroC,'xls','iDopaNeuroC')
        cd(currentDir)
    end
end
load([pathSave filesep 'iDopaNeuroC' filesep 'iDopaNeuroC' '.mat'])
eval(['model = ' 'iDopaNeuroC' ';'])
iDopaNeuroCStats = modelStatistics(model);
iDopaNeuroCData(7, 1) = iDopaNeuroCStats.Statistics(strcmp(iDopaNeuroCStats.Properties.RowNames,'Internal reactions'));
iDopaNeuroCData(8, 1) = iDopaNeuroCStats.Statistics(strcmp(iDopaNeuroCStats.Properties.RowNames,'External reactions'));
iDopaNeuroCData(9, 1) = iDopaNeuroCStats.Statistics(strcmp(iDopaNeuroCStats.Properties.RowNames,'Genes'));

%%

if 0
    %% compare two best thermoKernel models
    load([modelsDir filesep 'thermoKernel_oneRxnPerActiveGene_transcriptomicsT2_limitBoundary.10000_inactiveGenesT_closedIons_curationOverOmics' filesep 'Model.mat']);
    multimodels.model1=Model;
    load([modelsDir filesep 'thermoKernel_oneRxnPerActiveGene_transcriptomicsT2_limitBoundary.10000_inactiveGenesT_closedIons_omicsOverCuration' filesep 'Model.mat']);
    multimodels.model2=Model;
    [overlapresults,statistic] = compareXomicsModels(multimodels);
end

%% compare two best models
load([pathSave filesep 'iDopaNeuroC' filesep 'iDopaNeuroC'], 'iDopaNeuroC')
load([pathSave filesep 'iDopaNeuroCT' filesep 'iDopaNeuroCT'], 'iDopaNeuroCT')
multimodels.iDopaNeuroC=iDopaNeuroC;
multimodels.iDopaNeuroCT=iDopaNeuroCT;
[overlapresults,statistic] = compareXomicsModels(multimodels);
disp(statistic.overlaporportion_mets)
disp(statistic.overlaporportion_rxns)
disp(statistic.overlaporportion_genes)
disp(statistic.overlapnumber_mets)
disp(statistic.overlapnumber_rxns)
disp(statistic.overlapnumber_genes)

%% Print conditions
fid2 = fopen([pathSave filesep 'iDopaNeuroC' filesep 'conditions.txt'], 'w');
fprintf(fid2, '%s\n', iDopaNeuroCDir);
fclose(fid2);

%% Models statistics
% Table
bestModelsTable = table(iDopaNeuroCTData, iDopaNeuroCData, ...
    ...
    'VariableNames', ...
    {'iDopaNeuroCT', 'iDopaNeuroC'},...
    'RowNames',...
    {'Metabolites'; ...
    'Rank of S';...
    'Reactions';...
    'Qualitative';...
    'Spearman rho';...
    'Euclidean distance (μmol/gDW/hr)';...
    'Internal reactions';...
    'External reactions';...
    'Genes'});
display(bestModelsTable)

% Table figure
nexttile
% iDopaNeuroCT
ylocation = 0.1:(0.6/9):0.65;
text(0.04, ylocation(9), ' Cell-type (CT)', 'FontSize', 14, 'fontweight', 'bold');
text(0.18, ylocation(8), num2str(iDopaNeuroCTData(1)), 'FontSize', 14);
text(0.18, ylocation(7), num2str(round(iDopaNeuroCTData(3), 2)), 'FontSize', 14);
%text(0.18, ylocation(5), num2str(iDopaNeuroCTData(2)), 'FontSize', 14);
text(0.18, ylocation(6), num2str(round(iDopaNeuroCTData(7), 2)), 'FontSize', 14);
text(0.18, ylocation(5), num2str(round(iDopaNeuroCTData(8), 2)), 'FontSize', 14);
text(0.18, ylocation(4), num2str(round(iDopaNeuroCTData(9), 2)), 'FontSize', 14);
text(0.18, ylocation(3), num2str(round(iDopaNeuroCData(4), 2)), 'FontSize', 14);
text(0.18, ylocation(2), num2str(round(iDopaNeuroCTData(5), 2)), 'FontSize', 14);
text(0.18, ylocation(1), num2str(round(iDopaNeuroCTData(6), 2)), 'FontSize', 14);
% iDopaNeuroC
text(0.43, ylocation(9), ' Condition (C)', 'FontSize', 14, 'fontweight', 'bold');
text(0.43, ylocation(8), num2str(iDopaNeuroCData(1)), 'FontSize', 14);
text(0.43, ylocation(7), num2str(iDopaNeuroCData(3)), 'FontSize', 14);
text(0.43, ylocation(6), num2str(round(iDopaNeuroCData(7), 2)), 'FontSize', 14);
text(0.43, ylocation(5), num2str(round(iDopaNeuroCData(8), 2)), 'FontSize', 14);
text(0.43, ylocation(4), num2str(round(iDopaNeuroCData(9), 2)), 'FontSize', 14);
text(0.43, ylocation(3), num2str(round(iDopaNeuroCData(4), 2)), 'FontSize', 14);
text(0.43, ylocation(2), num2str(round(iDopaNeuroCData(5), 2)), 'FontSize', 14);
text(0.43, ylocation(1), num2str(round(iDopaNeuroCData(6), 2)), 'FontSize', 14);
yticks(ylocation)
set(gca,'YTickLabel', {'Euclidean distance', 'Spearman rho', 'Correct/total', 'Genes', 'Internal reactions', 'External reactions','rank(S)','Metabolites', 'iDopaNeuro'}, 'FontSize', 14, 'fontweight', 'bold')
title([letters(l) '. Statistics'], 'FontSize', 14)
l=l+1;
xline(0.4, 'LineWidth',3)
yline(0.599,'LineWidth',3)
axis([0 0.8 0.05 0.65])
set(gca,'XTick',[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose pairs of coloured labels in the plots
if 0
    idxConditionSpecific = conditionSpecificBool & modelsOverMinNumRxns & objectivesBool;
    idxCellTypeSpecific = cellTypeSpecificBool & modelsOverMinNumRxns & objectivesBool;
    legendLabels = {'Condition specific models', 'Cell-type specific models'};
elseif 1
    fastCoreBool = contains(multidimensionalComparisonStats.tissueSpecificSolver, 'fastCore');
    thermoKernelBool = contains(multidimensionalComparisonStats.tissueSpecificSolver, 'thermoKernel');
    idxConditionSpecific = fastCoreBool & modelsOverMinNumRxns & objectivesBool;
    idxCellTypeSpecific = thermoKernelBool & modelsOverMinNumRxns & objectivesBool;
    legendLabels = {'fastCore', 'thermoKernel'};
else
    unWeighted2normBool = ismember(multidimensionalComparisonStats.objective, 'unWeighted2norm');
    unWeightedTCBMfluxBool = ismember(multidimensionalComparisonStats.objective, 'unWeightedTCBMflux');
    idxConditionSpecific = unWeighted2normBool & modelsOverMinNumRxns & objectivesBool;
    idxCellTypeSpecific = unWeightedTCBMfluxBool & modelsOverMinNumRxns & objectivesBool;
    legendLabels = {'unWeighted2norm', 'unWeightedTCBMfluxBool'};
end

if allPlots
    %% Quantiative vs rank
    nexttile
    xFigF1 = multidimensionalComparisonStats.spearmanBoth(idxConditionSpecific);
    yFigF1 = multidimensionalComparisonStats.rankOfS(idxConditionSpecific);
    scatter(nonzeros(xFigF1), nonzeros(yFigF1), 'g')
    hold on
    xFigF2 = multidimensionalComparisonStats.spearmanBoth(idxCellTypeSpecific);
    yFigF2 = multidimensionalComparisonStats.rankOfS(idxCellTypeSpecific);
    scatter(nonzeros(xFigF2), nonzeros(yFigF2), 'mx')
    hold off
    title('D. Spearman rho vs rank(S)', 'FontSize', 14)
    xlabel('Spearman rho', 'fontweight', 'bold', 'FontSize', 14)
    ylabel('rank(S)', 'fontweight', 'bold', 'FontSize', 14)
    ymin = min([yFigF1; yFigF2], [], 'all');
    ymax = max([yFigF1; yFigF2], [], 'all');
    xmin = min([xFigF1; xFigF2],[],'all');
    xmax = max([xFigF1; xFigF2],[],'all');
    axis([xmin (xmax + (xmax * 0.05)) ymin ymax])
    text(multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCIdx), multidimensionalComparisonStats.rankOfS(iDopaNeuroCIdx), '• iDopaNeuroC', 'FontSize', 14);
    text(multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCTBool), multidimensionalComparisonStats.rankOfS(iDopaNeuroCTBool), '• iDopaNeuroCT', 'FontSize', 14);
end


%% Qualitative vs Spearman
nexttile
xFigE1 = multidimensionalComparisonStats.qualitativeBoth(idxConditionSpecific);
yFigE1 = multidimensionalComparisonStats.spearmanBoth(idxConditionSpecific);
scatter(nonzeros(xFigE1), nonzeros(yFigE1), 'g')
hold on
xFigE2 = multidimensionalComparisonStats.qualitativeBoth(idxCellTypeSpecific);
yFigE2 = multidimensionalComparisonStats.spearmanBoth(idxCellTypeSpecific);
scatter(nonzeros(xFigE2), nonzeros(yFigE2), 'mx')
hold off
title([letters(l) '. Correct/total vs Spearman'], 'FontSize', 14)
l=l+1;
xlabel('Correct/total predictions', 'fontweight', 'bold', 'FontSize', 14)
ylabel('Spearman rho', 'fontweight', 'bold', 'FontSize', 14)
ymin = min([yFigE1; yFigE2], [], 'all');
ymax = max([yFigE1; yFigE2], [], 'all');
xmin = min([xFigE1; xFigE2], [], 'all');
xmax = max([xFigE1; xFigE2], [], 'all');
axis([xmin (xmax + (xmax * 0.05)) ymin ymax])
text(multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCIdx), multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCIdx), '• iDopaNeuroC', 'FontSize', 14);
text(multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCTBool), multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCTBool), '• iDopaNeuroCT', 'FontSize', 14);
if legends
    legend(legendLabels, 'FontSize', 14)
end

if allPlots
    %% Qualitative vs Euclidean
    nexttile
    xFigE1 = multidimensionalComparisonStats.qualitativeBoth(idxConditionSpecific);
    yFigE1 = multidimensionalComparisonStats.quantitativeBoth(idxConditionSpecific);
    scatter(nonzeros(xFigE1), nonzeros(yFigE1), 'g')
    hold on
    xFigE2 = multidimensionalComparisonStats.qualitativeBoth(idxCellTypeSpecific);
    yFigE2 = multidimensionalComparisonStats.quantitativeBoth(idxCellTypeSpecific);
    scatter(nonzeros(xFigE2), nonzeros(yFigE2), 'mx')
    hold off
    title('F. Correct/total vs Euclidean', 'FontSize', 14)
    xlabel('Correct/total predictions', 'fontweight', 'bold', 'FontSize', 14)
    ylabel('Euclidean distance', 'fontweight', 'bold', 'FontSize', 14)
    ymin = min([yFigE1; yFigE2], [], 'all');
    ymax = max([yFigE1; yFigE2], [], 'all');
    xmin = min([xFigE1; xFigE2], [], 'all');
    xmax = max([xFigE1; xFigE2], [], 'all');
    axis([xmin (xmax + (xmax * 0.05)) ymin ymax])
    text(multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCIdx), multidimensionalComparisonStats.quantitativeBoth(iDopaNeuroCIdx), '• iDopaNeuroC', 'FontSize', 14);
    text(multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCTBool), multidimensionalComparisonStats.quantitativeBoth(iDopaNeuroCTBool), '• iDopaNeuroCT', 'FontSize', 14);
    if legends
        legend(legendLabels, 'FontSize', 14)
    end
end
set(gca,'fontsize',14,'FontWeight','bold')

%% Spearman vs Euclidean
nexttile
xFigE1 = multidimensionalComparisonStats.quantitativeBoth(idxConditionSpecific);
yFigE1 = multidimensionalComparisonStats.spearmanBoth(idxConditionSpecific);
scatter(nonzeros(xFigE1), nonzeros(yFigE1), 'g')
hold on
xFigE2 = multidimensionalComparisonStats.quantitativeBoth(idxCellTypeSpecific);
yFigE2 = multidimensionalComparisonStats.spearmanBoth(idxCellTypeSpecific);
scatter(nonzeros(xFigE2), nonzeros(yFigE2), 'mx')
hold off
title([letters(l) '. Euclidean vs Spearman'], 'FontSize', 14)
l=l+1;
xlabel('Euclidean distance', 'fontweight', 'bold', 'FontSize', 14)
ylabel('Spearman rho', 'fontweight', 'bold', 'FontSize', 14)
ymin = min([yFigE1; yFigE2], [], 'all');
ymax = max([yFigE1; yFigE2], [], 'all');
xmin = min([xFigE1; xFigE2], [], 'all');
xmax = max([xFigE1; xFigE2], [], 'all');
xmax = 50;%%%%%%%%%%%%%%%%%%%%%%%%%% just for paper
axis([xmin (xmax + (xmax * 0.05)) ymin ymax])
text(multidimensionalComparisonStats.quantitativeBoth(iDopaNeuroCIdx),multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCIdx), '• iDopaNeuroC', 'FontSize', 14);
text( multidimensionalComparisonStats.quantitativeBoth(iDopaNeuroCTBool),multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCTBool), '• iDopaNeuroCT', 'FontSize', 14);
if legends || 1
    legend(legendLabels, 'FontSize', 14)
end
set(gca,'fontsize',14,'FontWeight','bold')

if allPlots
    %% Qualitative vs rank
    nexttile
    xFigF1 = multidimensionalComparisonStats.qualitativeBoth(idxConditionSpecific);
    yFigF1 = multidimensionalComparisonStats.rankOfS(idxConditionSpecific);
    scatter(nonzeros(xFigF1), nonzeros(yFigF1), 'g')
    hold on
    xFigF2 = multidimensionalComparisonStats.qualitativeBoth(idxCellTypeSpecific);
    yFigF2 = multidimensionalComparisonStats.rankOfS(idxCellTypeSpecific);
    scatter(nonzeros(xFigF2), nonzeros(yFigF2), 'mx')
    hold off
    title('F. Correct/total vs rank(S)', 'FontSize', 14)
    xlabel('Correct/total predictions', 'fontweight', 'bold', 'FontSize', 14)
    ylabel('rank(S)', 'fontweight', 'bold', 'FontSize', 14)
    ymin = min([yFigF1; yFigF2], [], 'all');
    ymax = max([yFigF1; yFigF2], [], 'all');
    xmin = min([xFigF1; xFigF2],[],'all');
    xmax = max([xFigF1; xFigF2],[],'all');
    axis([xmin (xmax + (xmax * 0.05)) ymin ymax])
    text(multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCIdx), multidimensionalComparisonStats.rankOfS(iDopaNeuroCIdx), '• iDopaNeuroC', 'FontSize', 14);
    text(multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCTBool), multidimensionalComparisonStats.rankOfS(iDopaNeuroCTBool), '• iDopaNeuroCT', 'FontSize', 14);
end

%%
tmp = get(0, 'Screensize');
tmp(3)=tmp(3)/2;
set(gcf, 'Position', tmp);
savefig([pathSave filesep 'multidimensionalComparison'])
saveas(gcf,[pathSave filesep 'multidimensionalComparison'],'png')
saveas(gcf,[pathSave filesep 'multidimensionalComparison'],'eps')
saveas(gcf,[pathSave filesep 'multidimensionalComparison'],'jpg')