% Driver file to run AGORA2 pipeline on new organisms:
% Will run the six example draft reconstructions in
% semiautomaticrefinement/examples

% go to the folder where the draft reconstructions are
draftFolder = fileparts(which('example_infoFile.xlsx'));

% create a file with taxonomic information for your organism(s) similar to
% semiautomaticrefinement/examples/example_infoFile.xlsx and provide the path to the file
infoFilePath='example_infoFile.xlsx';

% set the input and output folders
mkdir('inputFiles')
inputDataFolder=[pwd filesep 'inputFiles' filesep];
addpath(inputDataFolder)

% define the name of the reconstruction project
reconVersion='Examples';
% if SBML files should be created
createSBML=false;

% if reports on each refined reconstruction in PDF format should be created.
% NOTE: Requires a working MikTex installation.
createReports=true;

% initialize the COBRA Toolbox and solvers
initCobraToolbox
% use IBM CPLEX if possible
solverOK=changeCobraSolver('ibm_cplex','LP');

% define number of workers for parallel computing
numWorkers=4;

% propagate input data if taxonomic information is provided
prepareInputData('infoFilePath', infoFilePath,'inputDataFolder', inputDataFolder)

%% 
% If species not yet included in AGORA2 are reconstructed: add data 
% manually at this point to the files in inputDataFolder based on avialable
% peer-reviewed publications or experimental data.
%%

% run pipeline
[refinedFolder, translDraftsFolder] = runPipeline(draftFolder, 'inputDataFolder', inputDataFolder, 'infoFilePath', infoFilePath, 'numWorkers', numWorkers, 'reconVersion', reconVersion, 'createSBML', createSBML);

% run test suite
runTestSuiteTools(translDraftsFolder, refinedFolder, 'numWorkers', numWorkers, 'infoFilePath', infoFilePath, 'inputDataFolder', inputDataFolder, 'reconVersion', reconVersion, 'createReports', createReports)

% Determine the distribution of properties across refined reconstructions
% if there are multiple ones (only if more than five strains are reconstructed)
computeModelProperties(translDraftsFolder, refinedFolder, 'numWorkers', numWorkers, 'infoFilePath', infoFilePath, 'reconVersion', reconVersion)
