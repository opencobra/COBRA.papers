%This script creates the variables through which the required parameters 
%and files are inputted to the metagenomic pipeline (MgPipe). Input 
%variables should be changed by the user according to what specified in the 
%documentation. Running this script will automatically launch the pipeline. 

% Federico Baldini, 2017-2018
% Almut Heinken, 08/2020: adapted to simplified inputs.

initCobraToolbox
global CBTDIR

%% REQUIRED INPUT VARIABLES

% path to microbe models
websave('AGORA-master.zip','https://github.com/VirtualMetabolicHuman/AGORA/archive/master.zip')
try
    unzip('AGORA-master')
end
modPath = [pwd filesep 'AGORA-master' filesep 'CurrentVersion' filesep 'AGORA_1_03' filesep' 'AGORA_1_03_mat'];

% path to and name of the file with abundance information.
abunFilePath=[CBTDIR filesep 'tutorials' filesep 'analysis' filesep 'microbiomeModelingToolbox' filesep 'normCoverageReduced.csv'];

% To define whether flux variability analysis to compute the metabolic profiles 
% should be performed
computeProfiles = true;

% path to and name of the file with dietary information
% (default='AverageEuropeanDiet')
dietFilePath = [CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'input' filesep 'AverageEuropeanDiet'];

% if to save models with diet constrains implemented (default=false)
saveConstrModels = true;

% number of cores dedicated for parallelization (default=2)
numWorkers = 4;

% Only inputs that you want to change from the default need to be declared.

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, computeProfiles, 'dietFilePath', dietFilePath, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers);

%% Statistical analysis and violin plots of the results
% Requires providing the path to a file with sample stratification
% information as the variable infoFilePath.
infoFilePath='sampInfo.csv';

% Header in the file with sample information with the stratification to 
% analyze (if not provided, the second column will be chosen by default)

sampleGroupHeaders={'Group'};
% sampleGroupHeaders can contain more than one entry if multiple columns 
% with sample information (e.g., disease state, age group) should be analyzed.

% to where the computation results to analyze are located
resPath = [pwd filesep 'Results'];

analyzeMgPipeResults(infoFilePath,resPath, 'sampleGroupHeaders', sampleGroupHeaders);

close all
