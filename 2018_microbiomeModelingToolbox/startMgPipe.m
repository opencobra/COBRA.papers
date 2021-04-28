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
system('curl -LJO https://github.com/VirtualMetabolicHuman/AGORA/archive/master.zip')
unzip('AGORA-master')
modPath = [pwd filesep 'AGORA-master' filesep 'CurrentVersion' filesep 'AGORA_1_03' filesep' 'AGORA_1_03_mat'];

% path to and name of the file with abundance information.
abunFilePath=[CBTDIR filesep 'tutorials' filesep 'analysis' filesep 'microbiomeModelingToolbox' filesep 'normCoverageReduced.csv'];

% To define whether flux variability analysis to compute the metabolic profiles 
% should be performed
computeProfiles = true;

% path where to save results (default=cobratoolbox/tmp)
mkdir('MicrobiomeModels')
resPath = [pwd filesep 'MicrobiomeModels'];

% path to and name of the file with dietary information
% (default='AverageEuropeanDiet')
dietFilePath = [CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'input' filesep 'AverageEuropeanDiet'];

% stratification of samples (note that group classification in the example 
% input file is not biologically meaningful)
infoFilePath='sampInfo.csv';

% to define whether the personalized mdoels should be pruned from a globals
% etup model (default), or created separately (recommended for large
% numbers of organisms)
buildSetupAll=true;

% if to save models with diet constrains implemented (default=false)
saveConstrModels = true;

% number of cores dedicated for parallelization (default=2)
numWorkers = 4;

% Only inputs that you want to change from the default need to be declared.

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'buildSetupAll', buildSetupAll, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers);

%% Statistical analysis and violin plots of the results
% Requires providing the path to a file with sample stratification
% information as the variable infoFilePath.
infoFilePath='sampInfo.csv';

% Header in the file with sample information with the stratification to 
% analyze (if not provided, the second column will be chosen by default)

sampleGroupHeaders={'Group'};
% sampleGroupHeaders can contain more than one entry if multiple columns 
% with sample information (e.g., disease state, age group) should be analyzed.

% define where results will be saved (optional, default folders will be
% generated otherwise)
statPath = [pwd filesep 'Statistics'];
violinPath = [pwd filesep 'ViolinPlots'];

analyzeMgPipeResults(infoFilePath,resPath,'statPath', statPath, 'violinPath', violinPath, 'sampleGroupHeaders', sampleGroupHeaders);

