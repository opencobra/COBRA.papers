% Generate iDopaNeuro models with different parameters using the
% function XomicsToMultipleModels. 
%
% The batch models are saved in 
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/multidimensionalModelGeneration


%% Set directory and specificData

clear
% specificData folder
dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics'];

%filenames
bibliomicData = 'bibliomicData.xlsx';
exometabolomicData = 'exometabolomicData.txt';
transcriptomicData = 'transcriptomicData.txt';

% generic model
inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' ...
    filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
genericModelName = 'Recon3DModel_301_xomics_input.mat';

%% automatic processing from here on
%bibliomic data
specificData = preprocessingOmicsModel([dataFolder filesep bibliomicData], 1, 1);

%exometabolomic data
specificData.exoMet = readtable([dataFolder filesep exometabolomicData]);

%transcriptomic data
specificData.transcriptomicData = readtable([dataFolder filesep transcriptomicData]);
specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);

replaceModels = false;

%load generic model
load([inputFolder filesep genericModelName])
modelGenerationConditions.genericModel.model = model;

if contains(char(java.lang.System.getProperty('user.name')),'rfleming')
    %action = 'debug';
    %action = 'oneModel';
    action = 'batch';
    action = 'iDopaNeuroC';
else
    action = 'batch';
end

switch action
    case 'batch'
        modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction'...
            filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'multidimensionalModelGeneration'];
        if ~isfolder(modelsDir)
            mkdir(modelsDir);
        end
%         cd(modelsDir)
        
        % Define results directory
        modelGenerationConditions.outputDir = modelsDir;
        % modelGenerationConditions
        % Note: if the data does not vary it is enough to declare them in param.
        % Here the non-varying conditions are left for demonstrative purposes
        %modelGenerationConditions.cobraSolver = {'gurobi'}; %  {'ibm_cplex', 'gurobi'};
        %modelGenerationConditions.cobraSolver = {'ibm_cplex'}; %  {'ibm_cplex', 'gurobi'};
        modelGenerationConditions.cobraSolver = {'mosek'}; 
        modelGenerationConditions.activeGenesApproach = {'deleteModelGenes', 'oneRxnPerActiveGene'}; % {'deleteModelGenes', 'oneRxnPerActiveGene'};
        modelGenerationConditions.transcriptomicThreshold = [0 2]; % [0 1 2];
        modelGenerationConditions.boundsToRelaxExoMet = {'b'}; % {'l', 'b', 'u'}
        modelGenerationConditions.closeIons = [true false]; % [true false];
        modelGenerationConditions.tissueSpecificSolver = {'fastCore', 'thermoKernel'}; % {'fastCore', 'thermoKernel'};
        modelGenerationConditions.inactiveGenesTranscriptomics = [true false]; % [true false];
        modelGenerationConditions.limitBounds = [1e4 1e5]; % [1e3 1e4 1e5];
        modelGenerationConditions.curationOverOmics = [true false]; % [true false];
        param.debug = false;
        param.printLevel = 1;
        
    case 'iDopaNeuroC'
        % thermoKernel_oneRxnPerActiveGene_transcriptomicsT2_limitBoundary.10000_inactiveGenesT_closedIons_omicsOverCuration

        modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
            filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'iDN1'];
        if ~isfolder(modelsDir)
            mkdir(modelsDir);
        end
        cd(modelsDir)
        
        % Define results directory
        modelGenerationConditions.outputDir = modelsDir;
        
        modelGenerationConditions.cobraSolver = {'ibm_cplex'};
        modelGenerationConditions.activeGenesApproach = {'oneRxnPerActiveGene'};
        modelGenerationConditions.transcriptomicThreshold = 2;
        modelGenerationConditions.boundsToRelaxExoMet = {'both'};
        modelGenerationConditions.closeIons = true;
        modelGenerationConditions.tissueSpecificSolver = {'thermoKernel'};
        modelGenerationConditions.inactiveGenesTranscriptomics = true; 
        modelGenerationConditions.limitBounds = 10000;
        modelGenerationConditions.curationOverOmics = false; 
        param.debug = true;
        param.printLevel = 2;
        
    case 'iDopaNeuroCT'
        % thermoKernel_oneRxnPerActiveGene_transcriptomicsT2_limitBoundary.10000_inactiveGenesT_closedIons_curationOverOmics

        modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
            filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'iDN1'];
        if ~isfolder(modelsDir)
            mkdir(modelsDir);
        end
        cd(modelsDir)
        
        % Define results directory
        modelGenerationConditions.outputDir = modelsDir;
        
        modelGenerationConditions.cobraSolver = {'ibm_cplex'};
        modelGenerationConditions.activeGenesApproach = {'oneRxnPerActiveGene'};
        modelGenerationConditions.transcriptomicThreshold = 2;
        modelGenerationConditions.boundsToRelaxExoMet = {'both'};
        modelGenerationConditions.closeIons = true;
        modelGenerationConditions.tissueSpecificSolver = {'thermoKernel'};
        modelGenerationConditions.inactiveGenesTranscriptomics = true; 
        modelGenerationConditions.limitBounds = 10000;
        modelGenerationConditions.curationOverOmics = true;
        param.debug = true;
        param.printLevel = 2;
    case 'debug'
        
        modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
            filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'iDopaNeuro1debug'];
        if ~isfolder(modelsDir)
            mkdir(modelsDir);
        end
        cd(modelsDir)
        
        % Define results directory
        modelGenerationConditions.outputDir = modelsDir;
        
        %debug generation of a single model type
        % modelGenerationConditions
        % Note: if the data does not vary it is enough to declare them in param.
        % Here the non-varying conditions are left for demonstrative purposes
        modelGenerationConditions.cobraSolver = {'gurobi'}; %  {'ibm_cplex', 'gurobi'};
        modelGenerationConditions.activeGenesApproach = {'deleteModelGenes'}; % {'deleteModelGenes', 'oneRxnPerActiveGene'};
        modelGenerationConditions.transcriptomicThreshold = 0; % [0 1 2];
        modelGenerationConditions.boundsToRelaxExoMet = {'both'}; % {'lower', 'both', 'upper'}
        modelGenerationConditions.closeIons = true; % [true false];
        modelGenerationConditions.tissueSpecificSolver = {'thermoKernel'}; % {'fastCore', 'thermoKernel'};
        modelGenerationConditions.inactiveGenesTranscriptomics = false; % [true false];
        modelGenerationConditions.limitBounds = 10000; % [1e3 1e4 1e5]; %was inf
        modelGenerationConditions.curationOverOmics = false; % [true false];
        param.debug = true;
        param.printLevel = 2;
end

% specificData
modelGenerationConditions.specificData.specificData = specificData;

%% Fixed options
param.setObjective = ''; % No objective function
feasTol = getCobraSolverParams('LP', 'feasTol');
param.boundPrecisionLimit = feasTol * 10;
param.fluxEpsilon = feasTol * 10;
param.fluxCCmethod = 'fastcc';
param.weightsFromOmics = 1;
param.metabolomicWeights = 'mean';
param.sinkDMinactive = 1; % Set non-core sinks and demands to inactive
param.addCoupledRxns = 1;
param.nonCoreSinksDemands = 'closeAll';
param.closeUptakes = true; % Cell culture information
param.findThermoConsistentFluxSubset = 1;


%% Create models

% Remove expressionRxns
if isfield(model, 'expressionRxns')
    model = rmfield(model, 'expressionRxns');
end

% Create models
directoriesWithModels = XomicsToMultipleModels(modelGenerationConditions, param);

disp('Models saved in:')
disp(modelsDir)
