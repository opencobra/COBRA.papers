% This file contains the minimal inputs needed to run the DEMETER pipeline.
% Adapt to your own organisms as necessary.

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

global CBTDIR

%% Preparing input data

% define the path to the folder where the draft reconstructions are located (required).
draftFolder = [CBTDIR filesep 'papers' filesep '2021_demeter' filesep 'exampleDraftReconstructions'];

% print names of refined reconstructions that will be created
refinedModelIDs = printRefinedModelIDs(draftFolder);

% define path to file with taxonomic information (create as needed)
infoFilePath = [CBTDIR filesep 'papers' filesep '2021_demeter' filesep 'example_infoFile.xlsx'];

% propagate experimental data
[infoFilePath,inputDataFolder] = prepareInputData(infoFilePath,'inputDataFolder', inputDataFolder);

%% Refining reconstructions

% Define the number of workers for parallel computing
numWorkers = 4;

% Define a name for the reconstruction resource (optional)
reconVersion = 'TutorialExample';

% Run the pipeline
[reconVersion,refinedFolder,translatedDraftsFolder,summaryFolder] = runPipeline(draftFolder, 'infoFilePath', infoFilePath, 'inputDataFolder', inputDataFolder,'numWorkers', numWorkers, 'reconVersion', reconVersion);

%% Testing reconstructions

testResultsFolder = runTestSuiteTools(refinedFolder, infoFilePath, inputDataFolder, reconVersion, 'translatedDraftsFolder', translatedDraftsFolder, 'numWorkers', numWorkers);

notGrowing = plotBiomassTestResults(refinedFolder,reconVersion,'translatedDraftsFolder',translatedDraftsFolder, 'numWorkers', numWorkers);

tooHighATP = plotATPTestResults(refinedFolder,reconVersion,'translatedDraftsFolder',translatedDraftsFolder, 'numWorkers', numWorkers);

%% Debugging the reconstructions

% Run the debugging suite
[debuggingFolder,debuggingReport, fixedModels, failedModels]=runDebuggingTools(refinedFolder,testResultsFolder,inputDataFolder,reconVersion,'numWorkers',numWorkers);

%% Plotting reconstructions

% Determine the distribution of properties across refined reconstructions.
% Will also cluster strains by similarity if enough different taxa are
% contained in the reconstruction resopurce.

propertiesFolder = computeModelProperties(refinedFolder, infoFilePath, reconVersion, 'numWorkers', numWorkers);

