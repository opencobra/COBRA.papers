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
system('curl -O https://www.vmh.life/files/reconstructions/AGORA/1.02/Agora-1.02.zip')
unzip('Agora-1.02.zip','AGORA')
modPath = [pwd filesep 'AGORA' filesep 'mat'];

% path to and name of the file with abundance information.
abunFilePath=[CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'examples' filesep 'normCoverage.csv'];

%% If you only want to set the required input variables, please run the
% pipeline as follows:
[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath);

% If you want to change any of the optional inputs, please find a
% description of them below.

%% OPTIONAL INPUTS
% path where to save results (default=cobratoolbox/tmp)
mkdir('MicrobiomeModels')
resPath = [pwd filesep 'MicrobiomeModels'];

% path to and name of the file with dietary information
% (default='AverageEuropeanDiet')
dietFilePath = [CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'resources' filesep 'AverageEuropeanDiet'];

% path to csv file for stratification criteria (if empty or not existent no criteria is used)
indInfoFilePath = '';

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
% flase=fluxVariability)
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

%% Pipeline start if setting any optional inputs
% Only inputs that you want to change from the default need to be declared.

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, 'resPath', resPath, 'dietFilePath', dietFilePath, 'indInfoFilePath', indInfoFilePath, 'objre', objre, 'figForm', figForm, 'numWorkers', numWorkers, 'autoFix', autoFix, 'compMod', compMod, 'rDiet', rDiet, 'pDiet', pDiet, 'extSolve', extSolve, 'fvaType', fvaType, 'autorun', autorun, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);


