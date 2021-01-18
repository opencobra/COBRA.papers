%This script creates the variables through which the required parameters 
%and files are inputted to the metagenomic pipeline (MgPipe). Input 
%variables should be changed by the user according to what specified in the 
%documentation. Running this script will automatically launch the pipeline. 

% Federico Baldini, 2017-2018
% Almut Heinken, 08/2020: adapted to simplified inputs.

initCobraToolbox()
global CBTDIR

%% REQUIRED INPUT VARIABLES

% path to microbe models
system('curl -LJO https://github.com/VirtualMetabolicHuman/AGORA/archive/master.zip')
unzip('AGORA-master')
modPath = [pwd filesep 'AGORA-master' filesep 'CurrentVersion' filesep 'AGORA_1_03' filesep' 'AGORA_1_03_mat'];

% path to and name of the file with abundance information.
abunFilePath=[CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'examples' filesep 'normCoverage.csv'];

%% If you only want to set the required input variables, please run the
% pipeline as follows:
[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath);

% If you want to change any of the optional inputs, please find a
% description of them below.

%% Pipeline start if setting any optional inputs

% path where to save results (default=cobratoolbox/tmp)
mkdir('MicrobiomeModels')
resPath = [pwd filesep 'MicrobiomeModels'];

% path to and name of the file with dietary information
% (default='AverageEuropeanDiet')
dietFilePath = [CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'resources' filesep 'AverageEuropeanDiet'];

% stratification of samples (note that group classification in the example 
% input file is not biologically meaningful)
infoFilePath='sampInfo.csv';

% name of objective function of organisms, default='EX_biomass(e)'
objre = 'EX_biomass(e)';

%the output is vectorized picture, default=-depsc, change to '-dpng' for .png
figForm = '-depsc';

% number of cores dedicated for parallelization (default=2)
numWorkers = 2;

% autofix for names mismatch (default=true)
autoFix = true;

% if outputs in open formats should be produced for each section (default=false)
compMod = false; 

% to enable also rich diet simulations (default=false)
rDiet = false;

% to enable personalized diet simulations (default=false)
pDiet = false;

% if to use an external solver and save models with diet (default=false)
extSolve = false;

% the type of FVA function to use to solve (true=fastFVA,
% false=fluxVariability)
fvaType = true;

% To turn off the autorun to be able to manually execute each part of the pipeline (default=true).
autorun = true;

% to manually set the lower bound on flux through the community biomass
% reaction (default=0.4 mmol/person/day)
lowerBMBound = 0.4;

% to set whether existing simulation results are rewritten (default=false)
repeatSim = false;

% to set if the input medium should be adapted through the adaptVMHDietToAGORA
% function or used as is (default=true)                  
adaptMedium = true; 

% Only inputs that you want to change from the default need to be declared.

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'objre', objre, 'figForm', figForm, 'numWorkers', numWorkers, 'autoFix', autoFix, 'compMod', compMod, 'rDiet', rDiet, 'pDiet', pDiet, 'extSolve', extSolve, 'fvaType', fvaType, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);

%% Pipeline start if including Recon3D as the host

system('curl -LJO https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D_301.zip')
unzip('Recon3D_301')
hostPath = [pwd filesep 'Recon3D_301' filesep 'Recon3DModel_301.mat'];
%
% Since host metabolites can now enter from the host model itself, the 
% adaptMedium input can be set to false.                 
adaptMedium = false; 
%
% overwrite the previous simulation for the purpose of the tutorial
repeatSim = true;
%
% If a host model is entered, it is also highly recommended to enter the
% host biomass reaction to generate coupling constraints for the host.
hostBiomassRxn = 'biomass_reaction';

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'hostPath', hostPath, 'hostBiomassRxn', hostBiomassRxn, 'objre', objre, 'figForm', figForm, 'numWorkers', numWorkers, 'autoFix', autoFix, 'compMod', compMod, 'rDiet', rDiet, 'pDiet', pDiet, 'extSolve', extSolve, 'fvaType', fvaType, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);

%% Statistical analysis and violin plots of the results
% Requires providing the path to a file with sample stratification
% information as the variable infoFilePath.
infoFilePath='sampInfo.csv';

% Header in the file with sample information with the stratification to 
% analyze (if not provided, the second column will be chosen by default)

sampleGroupHeaders={'Group'};
% sampleGroupHeaders can contain more than one entry if multiple columns 
% with sample information (e.g., disease state, age group) should be analyzed.

% path with results of mgPipe that will be analyzed
resPath = [tutorialPath filesep 'Results'];

% define where results will be saved (optional, default folders will be
% generated otherwise)
statPath = [tutorialPath filesep 'Statistics'];
violinPath = [tutorialPath filesep 'ViolinPlots'];

analyzeMgPipeResults(infoFilePath,resPath,'statPath', statPath, 'violinPath', violinPath, 'sampleGroupHeaders', sampleGroupHeaders);

