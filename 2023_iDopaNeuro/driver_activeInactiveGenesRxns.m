% Generates a graph that shows the bibliomics active genes and reactions in
% the iDopaNeuro models
%
% The graph is saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1/

clear

% specificData
dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];
bibliomicData = 'bibliomicData.xlsx';

% Generic model
inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' ...
    filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
genericModelName = 'Recon3DModel_301_xomics_input.mat';
load([inputFolder filesep genericModelName])
genericModel = model;

% Select a COBRA model
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

if 0
    models = {'iDopaNeuroCT'; 'iDopaNeuroC'};
else
    models = {'iDopaNeuroC'};
end

%% Genes

    tableGenes = readtable([dataFolder bibliomicData], 'Sheet', 'activeGenes');
    tableGenes.source(contains(tableGenes.source, 'Dopaminergic neuron manual curation (metabolism)')) = {'Dopamine metabolism'};
    tableGenes.source(contains(tableGenes.source, 'Dopaminergic neuron manual curation (transport genes)')) = {'Cellular transport'};
    tableGenes.source(contains(tableGenes.source, 'mitochondrial manual curation from (Elassal.D et al. In preparation)')) = {'Mitochondrial metabolism'};
    tableGenes.source(contains(tableGenes.source, 'manual curation from (Elassal.D et al. In preparation)')) = {'Carbon central metabolism'};
    
    % Remove bibliomic genes not present in the model
    tableGenes(~ismember(tableGenes.geneID, str2double(regexprep(genericModel.genes, '\.\d', ''))), :) = [];
    
%% Reactions
    
    % Read activeReactions
    tableRxns = readtable([dataFolder bibliomicData], 'Sheet', 'activeReactions');
    tableRxns = removevars(tableRxns, setdiff(tableRxns.Properties.VariableNames, {'rxns', 'type'}));
    % Read sinkDemand
    % tableRxns2 = readtable([dataFolder bibliomicData], 'Sheet', 'sinkDemand');
    % tableRxns2 = removevars(tableRxns2, setdiff(tableRxns2.Properties.VariableNames, {'rxns'}));
    % tableRxns2.type(1:size(tableRxns2, 1)) = {'Sink and demands'};
    % [~, ic, ia] = unique(tableRxns2.rxns);
    % tableRxns2(setdiff(ia, ic), :) = [];
    % Read Medium mets
    tableRxns3 = readtable([dataFolder bibliomicData], 'Sheet', 'mediaData');
    tableRxns3 = removevars(tableRxns3, setdiff(tableRxns3.Properties.VariableNames, {'rxns'}));
    tableRxns3.type(1:size(tableRxns3, 1)) = {'Media metabolites'};
    tableRxns3(ismember(tableRxns3.rxns, tableRxns.rxns), :) = [];
    % Read Biomass precursors
    tableRxns4 = readtable([dataFolder bibliomicData], 'Sheet', 'coupledRxns');
    rxns = [];
    for i = 1:length(tableRxns4.coupledRxnsList)
        rxns = [rxns; split(tableRxns4.coupledRxnsList{i}, ', ')];
    end
    tableRxns4 = table(rxns);
    tableRxns4.type(1:size(tableRxns4, 1)) = {'Biomass precursors'};
    tableRxns(ismember(tableRxns.rxns, tableRxns4.rxns), :) = [];
    % Join tables
    % tableRxns = [tableRxns; tableRxns2; tableRxns3; tableRxns4];
    tableRxns = [tableRxns; tableRxns3; tableRxns4];
    % Read rxns2constrain
    tableRxns5 = readtable([dataFolder bibliomicData], 'Sheet', 'rxns2constrain');
    rxns = tableRxns5.rxns;
    tableRxns5 = table(rxns);
    tableRxns5.type(1:size(tableRxns5, 1)) = {'Others'};
    tableRxns5(ismember(tableRxns5.rxns, tableRxns.rxns), :) = [];
    % Join tables
    tableRxns = [tableRxns; tableRxns5];
    tableRxns.type(cellfun(@isempty, tableRxns.type)) = {'Others'};
    
    % Rename type
    tableRxns.type(contains(tableRxns.type, 'Other')) = {'Others'};
    tableRxns.type(contains(tableRxns.type, 'TranspMets+DiffMediumMets')) = {'Cellular transport'};
    tableRxns.type(contains(tableRxns.type, 'Dopaminergic neuron metabolism')) = {'Dopamine metabolism'};
    tableRxns.type(contains(tableRxns.type, 'manual curation from (Elassal.D et al. In preparation)')) = {'Mitochondrial & carbon central metabolism'};
    tableRxns.type(contains(tableRxns.type, 'Sinks and demands')) = {'Others'};
    tableRxns.type(contains(tableRxns.type, 'Exometabolomics [completely secreted metabolites – not ini')) = {'Others'};
    
    % Remove bibliomic genes not present in the model
    tableRxns(~ismember(tableRxns.rxns, genericModel.rxns), :) = [];
    
%% Pie Charts

for j = 1:length(models)
    
    % Load model
    load([pathSave filesep models{j} filesep models{j} '.mat'])
    eval(['model = ' models{j} ';'])
    
    % GENES
    
    % Figure
    % Count labels
    tableGenesExcluded = tableGenes;
    modelGenes = str2double(regexprep(model.genes, '\.\d', ''));
    activeGenes = tableGenesExcluded.geneID;
    tableGenesExcluded.excluded = ~ismember(activeGenes, modelGenes);
    tableGenesExcluded.source(~ismember(activeGenes, modelGenes)) = ...
        strcat(tableGenes.source(~ismember(activeGenes, modelGenes)), ' (excluded)');
    [uniqueLabelsGenes, ~, idxGenes] = unique(tableGenesExcluded.source, 'stable');
    count = hist(idxGenes, unique(idxGenes));
    % Split excluded genes
    [count, ia] = sort(count);
    uniqueLabelsGenes = uniqueLabelsGenes(ia);
    excludedIdx = contains(uniqueLabelsGenes, 'excluded');
    countExcluded = count(excludedIdx);
    count = [count(~excludedIdx) sum(countExcluded)];
    uniqueLabelsGenesExcluded = uniqueLabelsGenes(excludedIdx);
    uniqueLabelsGenes = [uniqueLabelsGenes(~excludedIdx); {'Excluded'}];
    
    % Pie chart
    figure
    subplot(2, 2, 1)
    pieChart = pie(count);
    ax = gca();
    newColors = [...
        0         0.4471    0.7412;
        0.8510    0.3255    0.0980;
        0.9294    0.6941    0.1255;
        0.4941    0.1843    0.5569;
        0.1490    0.1490    0.1490];
    ax.Colormap = newColors;
    title(['A. Bibliomics active genes: ' num2str(sum(count))], 'FontSize', 14)
    legend(uniqueLabelsGenes, 'FontSize', 14, 'Location', 'best');
    
    % Excluded
    subplot(2, 2, 2)
    pieChart = pie(countExcluded);
    ax = gca();
    rows = ismember(uniqueLabelsGenes, regexprep(uniqueLabelsGenesExcluded, ...
        ' \(excluded\)', ''));
    newColors = newColors(rows, :);
    ax.Colormap = newColors;
    title(['B. Excluded genes: ' num2str(sum(countExcluded))], 'FontSize', 14)
    % legend(uniqueLabelsGenes(rows), 'FontSize', 12, 'Location', 'best');
    
    tableGenesExcluded = tableGenesExcluded(tableGenesExcluded.excluded==1,1);
    
    % REACTIONS
    
    % Figure
    % Count labels
    tableRxnsExcluded = tableRxns;
    modelRxns = model.rxns;
    activeRxns = tableRxnsExcluded.rxns;
    tableRxnsExcluded.excluded = ~ismember(activeRxns, modelRxns);
    tableRxnsExcluded.type(~ismember(activeRxns, modelRxns)) = ...
        strcat(tableRxns.type(~ismember(activeRxns, modelRxns)), ' (excluded)');
    [uniqueLabelsRxns, ~, idxRxns] = unique(tableRxnsExcluded.type, 'stable');
    count = hist(idxRxns, unique(idxRxns));
    % Split excluded genes
    [count, ia] = sort(count);
    uniqueLabelsRxns = uniqueLabelsRxns(ia);
    excludedIdx = contains(uniqueLabelsRxns, 'excluded');
    countExcluded = count(excludedIdx);
    count = [count(~excludedIdx) sum(countExcluded)];
    uniqueLabelsRxnsExcluded = uniqueLabelsRxns(excludedIdx);
    uniqueLabelsRxns = [uniqueLabelsRxns(~excludedIdx); {'Excluded'}];
    
    % Pie chart
    subplot(2, 2, 3)
    pieChart = pie(count);
    ax = gca();
    newColors = [...
        0.6353    0.0784    0.1843; % MM brown
        0.8510    0.3255    0.0980; % CT red
        0.7176    0.2745    1.0000; % BP purple
        1.0000    0.0745    0.6510; % OT pink
        0         0.4471    0.7412; % DM blue
        0.9294    0.6941    0.1255; % CCM yellow
        0.1490    0.1490    0.1490];% EXC black
    ax.Colormap = newColors;
    title(['C. Bibliomics active reactions: ' num2str(sum(count))], 'FontSize', 14)
    legend(uniqueLabelsRxns, 'FontSize', 14, 'Location', 'best');
    
    % Excluded
    subplot(2, 2, 4)
    [bool, locb] = ismember(regexprep(uniqueLabelsRxnsExcluded,' \(excluded\)', ''), uniqueLabelsRxns);
    pieData = zeros(size(count));
    pieData(locb) = countExcluded;
    pieChart = pie(pieData);
    ax = gca();
    ax.Colormap = newColors;
    title(['D. Excluded reactions: ' num2str(sum(countExcluded))], 'FontSize', 14)
    % legend(uniqueLabelsRxns(rows), 'FontSize', 12, 'Location', 'best');
    
    tableRxnsExcluded = tableRxnsExcluded(tableRxnsExcluded.excluded==1,1);
    % Save figure
    savefig([pathSave filesep models{j} filesep 'bibliomicsGenesRxns_' models{j}])
    saveas(gcf,[pathSave filesep models{j} filesep 'bibliomicsGenesRxns_' models{j}],'png')
    saveas(gcf,[pathSave filesep models{j} filesep 'bibliomicsGenesRxns_' models{j}],'eps')
    
    clearvars -except j models tableRxns tableGenes pathSave tableRxnsExcluded tableGenesExcluded
end