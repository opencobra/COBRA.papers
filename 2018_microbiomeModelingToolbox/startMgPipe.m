%This script creates the variables through which the required parameters 
%and files are inputted to the metagenomic pipeline (MgPipe). Input 
%variables should be changed by the user according to what specified in the 
%documentation. Running this script will automatically launch the pipeline. 

% Federico Baldini, 2017-2018

initCobraToolbox()

%% REQUIRED INPUT VARIABLES

% path to microbe models
system('curl -LJO https://github.com/VirtualMetabolicHuman/AGORA/archive/master.zip')
unzip('AGORA-master')
modPath = [pwd filesep 'AGORA-master' filesep 'CurrentVersion' filesep 'AGORA_1_03' filesep' 'AGORA_1_03_mat'];

% path to and name of the file with abundance information.
abunFilePath='normCoverage.csv';

%% OPTIONAL INPUT VARIABLES

% path where to save results
mkdir('MicrobiomeModels')
resPath=[pwd filesep 'MicrobiomeModels'];

% path to and name of the file with dietary information.
dietFilePath='AverageEuropeanDiet';

% number of cores dedicated for parallelization 
numWorkers = 3;


%% PIPELINE LAUNCHER 
[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, 'resPath', resPath, 'dietFilePath', dietFilePath, 'numWorkers', numWorkers);

