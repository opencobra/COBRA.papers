% Specific data preprocessing
% The tables transcriptomicsData and exoMet, which are used in the variable
% specificData, are generated using the function transcriptomics2genes 

% The external files are saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/data/xomics

clear

dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];

%% transcriptomicData
%   * .genes
%   * .expVal

if ~isfile([dataFolder 'transcriptomicData.txt'])
    
    display('Processing transcriptomic data ...')
    transcriptomicFile = 'SM7.tab';
    transcriptomicData = transcriptomics2genes([dataFolder transcriptomicFile]);
    % Save table
    writetable(transcriptomicData, [dataFolder 'transcriptomicData'])
    display('The transcriptomic data is ready')
    
else
    display('The transcriptomic data is ready')
end

%% exoMet
%   * .rxns
%   * .rxnNames
%   * .mean
%   * .SD
%   * .units
%   * .Platform

if ~isfile([dataFolder 'exometabolomicData.txt'])
    
    display('Processing exometabolomic data ...')
    
    exoMetData = '~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/metabolomicAnalysis/combinedRates.mat';
    load(exoMetData)
    
    % Prepare the table with the exometabolomic data
    nRows = size(combinedData, 1);
    varTypes = {'string', 'string', 'double', 'double', 'string', 'string'};
    varNames = {'rxns', 'rxnNames', 'mean', 'SD', 'units', 'platform'};
    exoMet = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes, 'VariableNames', varNames);
    
    % Create table based on Cornelius data
    exoMet.rxns = combinedData.rxns;
    exoMet.rxnNames = append('Exchange of ', combinedData.metName);
    exoMet.mean = combinedData.mean_ctrl;
    exoMet.SD = combinedData.sd_ctrl;
    exoMet.units(:) = {'μmol/gDW/hr'};
    exoMet.platform = combinedData.platform;
    % Delete NaN data
    exoMet(isnan(exoMet.SD), :) = [];
    
%     % Create table based on Can data
%     conc2Fluxes
%     rxns = v(1, 2:end)';
%     rxnNames = cell(size(rxns));
%     mean = [v{11, 2:end}]';
%     SD = [v{12, 2:end}]';
%     units = cell(size(rxns));
%     units(:) = {'μmol/gDW/hr'};
%     Platform = v(17, 2:end)';
%     exoMet2 = table(rxns, rxnNames, mean, SD, units, Platform);
%     % Preference to Cornelius data
%     repeatedMetsBool = ismember(exoMet2.rxns, exoMet.rxns);
%     exoMet2(repeatedMetsBool, :) = [];
    
%     % Merge and save
%     exoMet = [exoMet; exoMet2];
    
    % Save table
    dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
        filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics'];
    writetable(exoMet, [dataFolder filesep 'exoMet'])
    display(['The exometabolomic data is saved as:' dataFolder filesep 'exometabolomicData.txt'])
    
else
    display(['The exometabolomic data already exists as:' dataFolder filesep 'exometabolomicData.txt'])
end