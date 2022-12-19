% Generates a venn diagram that compares the core reactions between the
% iDopaNeuroCT and iDopaNeuroC models. Additionally, it is compared the 
% information of both models such as dimentions, essential pathways, main
% energy sources, etc.
%
% The results are saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1

%% Define directories

clear

% Results dir
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

% specificData dir
dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];
bibliomicData = 'bibliomicData.xlsx';

%% Core reactions

% Read activeReactions
specificData = preprocessingOmicsModel([dataFolder bibliomicData], 1, 1);
% coreRxnAbbr
coreRxnAbbr = {};
coreRxnAbbr = [coreRxnAbbr; specificData.activeReactions];
coreRxnAbbr = [coreRxnAbbr; specificData.rxns2constrain.rxns];
for i = 1:length(specificData.coupledRxns.coupledRxnsList)
    coreRxnAbbr = [coreRxnAbbr; split(specificData.coupledRxns.coupledRxnsList{i}, ', ')];
end
coreRxnAbbr = [coreRxnAbbr; specificData.mediaData.rxns];
coreRxnAbbr = unique(coreRxnAbbr);

%% Venn diagram

% Load models
load([pathSave filesep 'iDopaNeuroCT' filesep 'iDopaNeuroCT.mat'])
load([pathSave filesep 'iDopaNeuroC' filesep 'iDopaNeuroC.mat'])

iDopaNeuroCT_rxns = iDopaNeuroCT.rxns;
iDopaNeuroC_rxns = iDopaNeuroC.rxns;

%% Generate diagram
set7 = numel(intersect(iDopaNeuroC_rxns, intersect(iDopaNeuroCT_rxns, coreRxnAbbr)));
vennX([numel(setdiff(coreRxnAbbr, [iDopaNeuroCT_rxns; iDopaNeuroC_rxns])) ...
    numel(intersect(coreRxnAbbr, iDopaNeuroCT_rxns))  - set7 ...
    numel(setdiff(iDopaNeuroCT_rxns, [coreRxnAbbr; iDopaNeuroC_rxns]))  ...
    numel(intersect(iDopaNeuroC_rxns, iDopaNeuroCT_rxns))  - set7 ...
    numel(setdiff(iDopaNeuroC_rxns, [iDopaNeuroCT_rxns; coreRxnAbbr]))  ...
    numel(intersect(coreRxnAbbr, iDopaNeuroC_rxns))  - set7 ...
    set7], .01);
title({'Core reactions between', 'iDopaNeuro models'}, 'FontSize', 16)

% Save figure
savefig([pathSave filesep 'coreRxnsComparisonVenn'])
saveas(gcf,[pathSave filesep 'coreRxnsComparisonVenn'],'png')
%saveas(gcf,[pathSave filesep 'coreRxnsComparisonVenn'],'eps')

%% iDopaNeuro models comparison

% Select the models to be compared
rxnsList = {'EX_atp[e]', 'EX_dopa[e]', 'EX_dopasf[e]', 'EX_h2o2[e]', 'EX_dopa4sf[e]', ...
    'EX_dopa4glcur[e]', 'EX_dopa3glcur[e]', 'EX_4glu56dihdind[e]', 'EX_5cysdopa[e]', 'EX_CE5025[e]'};
iDopaNeuroC = changeRxnBounds(iDopaNeuroC, rxnsList, 0, 'l');
iDopaNeuroCT = changeRxnBounds(iDopaNeuroCT, rxnsList, 0, 'l');

models.iDopaNeuroC = iDopaNeuroC;
models.iDopaNeuroCT = iDopaNeuroCT;

% Select the type of comparison preformed
param.objectives = {'unWeightedTCBMflux'};
param.trainingSet = iDopaNeuroCT.XomicsToModelSpecificData.exoMet;
param.printLevel = 0;
param.tests = {'flux'};

[solverOK, solverInstalled] = changeCobraSolver('mosek','LP', 0);

%% Robustness Analysis

% Define parameters
objective = {'unWeightedTCBMflux'};
param.objectives = objective;
param.printLevel = 0;
param.tests = {'flux'};
rxn2vary = 'ATPM';
fixedRxns = {'PGK'; 'PYK'; 'ATPtm'; 'GLUDxm'};
range = 186:186 * 10;
npoints = 100;
bound = 'l';

% find fixedRxns idx
iDNc_rxn2Idxs = findRxnIDs(iDopaNeuroC, fixedRxns);
iDNct_rxn2Idxs = findRxnIDs(iDopaNeuroCT, fixedRxns);

% Compute fluxes
[rx1, rx2] = deal(zeros(npoints, 1));
bounds = min(range):(max(range) - min(range)) / npoints:max(range);
for i = 1:length(fixedRxns)
    
    for j = 1:length(bounds)
        
        % iDopaNeuroC
        iDNc = iDopaNeuroC;
        iDNc = changeRxnBounds(iDNc, rxn2vary, bounds(j), bound);
        solutions = modelMultipleObjectives(iDNc, param);
        iDNc_rx1(j) = bounds(j);
        iDNc_rx2(j, i) = solutions.(objective{:}).v(iDNc_rxn2Idxs(i));
        
        % iDopaNeuroC
        iDNct = iDopaNeuroCT;
        iDNct = changeRxnBounds(iDNct, rxn2vary, bounds(j), bound);
        solutions = modelMultipleObjectives(iDNct, param);
        iDNct_rx1(j) = bounds(j);
        iDNct_rx2(j, i) = solutions.(objective{:}).v(iDNct_rxn2Idxs(i));
        
    end
end
save([pathSave filesep 'iDopaNeuroModelsComparison.mat'])

%% Generate figures
for i = 1:length(fixedRxns)
    switch i
        case 1
            label = 'B.';
        case 2
            label = 'C.';
        case 3
            label = 'D.';
        case 4
            label = 'E.';
    end
    sgtitle({'Phenotype of essential reactions when varying ATP demand'}, 'FontSize', 16)
    subplot(2, 2, i)
    plot(iDNc_rx1, iDNc_rx2(:, i), 'LineWidth', 2)
    hold on 
    plot(iDNct_rx1, iDNct_rx2(:, i), 'LineWidth', 2, 'Color', [0.9294 0.6941 0.1255])
    axis([186 186*10 min(min([iDNc_rx2(:, i) iDNct_rx2(:, i)])) max(max([iDNc_rx2(:, i) iDNct_rx2(:, i)]))])
    xlabel([rxn2vary ' umol/gDW/hr'], 'FontSize', 12)
    ylabel([fixedRxns{i} ' umol/gDW/hr'], 'FontSize', 12,'Interpreter', 'none')
    title([label ' ' iDopaNeuroCT.rxnNames{iDNct_rxn2Idxs(i)}], 'FontSize', 14)
    if i == 1
        legend({'iDopaNeuroC'; 'iDopaNeuroCT'}, 'Location', 'best', 'FontSize', 12)
    end
end
%%

savefig([pathSave filesep 'robustnessAnalysis'])
saveas(gcf,[pathSave filesep 'robustnessAnalysis'],'png')
%saveas(gcf,[pathSave filesep 'rba'],'eps')

