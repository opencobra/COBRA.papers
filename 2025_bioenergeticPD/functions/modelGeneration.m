%clear all
clear
beep on
% initCobraToolbox
%%
% preconditional recon model
% preConditionRecon_bioenergetics
%% 1. genericModel
% genericModelName = 'Recon3D_301.mat'; 
% genericModelName = 'Recon3DModel_301.mat'; % recon3DModel
genericModelName = 'Recon3DModel_301_xomics_input.mat'; % using the xomics input
% genericModelName = 'Recon3d.mat' ;% Preconditional recon model
% genericModelName = 'Recon3dnew.mat' ;% with revised FAD-rxns
% genericModelName = 'Recon3dunlumped.mat' ;% with revised FAD-rxns + unlumped rxns
switch genericModelName
    case 'Recon3DModel_301.mat'
        load('~/fork-cobratoolbox/test/models/mat/Recon3DModel_301.mat');
        modelOrig=model;
    case 'Recon3D_301.mat'
         load('~/fork-cobratoolbox/test/models/mat/Recon3D_301.mat');
        model=Recon3D;
    case 'Recon3DModel_301_xomics_input.mat'
        load('~/drive/sbgCloud/projects/xomicsToModel/data/Recon3D_301/Recon3DModel_301_xomics_input.mat');
        modelOrig=model;
    case 'Recon3d.mat'
        %load('~/drive/bioenergeticsPD/rawData/Recon3d.mat')
        load('~/drive/bioenergeticsPD/fromXi/Recon3d.mat')
        modelOrig=model;
    case 'Recon3dnew.mat'
        load('~/drive/bioenergeticsPD/fromXi/Recon3dnew.mat')
        modelOrig=model;
    case 'Recon3dunlumped.mat'
         load('~/drive/bioenergeticsPD/fromXi/Recon3dunlumped.mat')
        modelOrig=model;
end
%% 1.1 revise GBA gene
GPR=model.grRules{ismember(model.rxns,'GBA2e')};% GBA1 gene only in lysosome, so change the cytosol gene to GBA2
model = changeGeneAssociation(model, 'GBA', GPR);

%% 1.2 refine the FAD rxns on Recon3D (skip)
% % add new FAD rxns
% % FAD=readtable('~/drive/bioenergeticsPD/fromXi/data/FDHrxns/FAD_Correction.xlsx','ReadVariableNames',true);
%% 1.3 select threshold for the transcriptomoc data
SelectThreshold
%% 2. Context-specific data
dataFolder = '~/drive/bioenergeticsPD/fromXi/data/inputData';

%%%% add GEO and nonGEO results in the reconstruction file
% bibliomicData = 'synaptic_bibliomicData.xlsx';
% bibliomicData = 'synapticPD_bibliomicData.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData_unconstrained.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData_unconstrained.xlsx';
% bibliomicData = 'astro_bibliomicData.xlsx';
% bibliomicData = 'astro_bibliomicData_unconstrained.xlsx';
% specificData = preprocessingOmicsModel([dataFolder filesep bibliomicData], 1, 1);

%%%% new3: only raw reconstruction file from Diana (without GEO and nonGEO active genes)
% bibliomicData = 'synaptic_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained_reconstruction.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained_reconstruction.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData_unconstrained_reconstruction.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData_unconstrained_reconstruction.xlsx';
% specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep 'new3' filesep bibliomicData], 1, 1);

%%%% newly organized biblomics data using exometabolomic data from iDopa
% bibliomicData = 'test_SYN.xlsx'; % revised rxns2constrain file for synaptic_bibliomicData_reconstruction.xlsx
% bibliomicData = 'test_SYNPD.xlsx';
% bibliomicData = 'test_ASYN.xlsx';
% bibliomicData = 'test_ASYNPD.xlsx';

% bibliomicData = 'test_SYN_new.xlsx';
bibliomicData = 'test_SYNPD_new.xlsx'; 
% bibliomicData = 'test_ASYN_new.xlsx';
% bibliomicData = 'test_ASYNPD_new.xlsx';
specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep 'new3' filesep bibliomicData], 1, 1);

%%%% raw reconstruction file from Diana (without GEO and nonGEO active
%%%% genes) + newly added active rxns
% bibliomicData = 'synaptic_bibliomicData_reconstruction1.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_reconstruction1.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained_reconstruction1.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained_reconstruction1.xlsx';

% specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep 'new4' filesep bibliomicData], 1, 1);
%% use GEO transcriptomic expression value
% dataFile='~/drive/bioenergeticsPD/fromXi/data/snRNA-seq/GSE157783/data';
dataFile='~/drive/bioenergeticsPD/fromXi/data/snRNA-seq/GSE178265/CellXgene';
if any(strfind(bibliomicData,'synapticPD')) & ~any(strfind(bibliomicData,'astro'))
%     specificData.transcriptomicData = readtable([dataFolder filesep 'GEO_transcriptomicData_allPD_geneval.txt']);
%     specificData.transcriptomicData = readtable([dataFolder filesep 'GEO_transcriptomicData_allPD_zscore.txt']);
%     specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);
% snRNA-seq
%     specificData.transcriptomicData = readtable([dataFile filesep 'SCTnorm_pd.csv']);
%       specificData.transcriptomicData = readtable([dataFile filesep 'scaled_pd.csv']);
%     specificData.transcriptomicData = readtable([dataFile filesep 'DNpd_norexp1.csv']);
    % specificData.transcriptomicData = readtable([dataFile filesep 'DNpd_NB.csv']);
    specificData.transcriptomicData = readtable([dataFile filesep 'DNpd_NBnew.csv']);
    param.transcriptomicThreshold=-3; % for both DNpd_NB.csv and DNpd_NBnew.csv
    
elseif ~any(strfind(bibliomicData,'synapticPD')) & ~any(strfind(bibliomicData,'astro'))
%     specificData.transcriptomicData = readtable([dataFolder filesep 'GEO_transcriptomicData_allControl_geneval.txt']);
%     specificData.transcriptomicData = readtable([dataFolder filesep 'GEO_transcriptomicData_allControl_zscore.txt']);
%     specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);
% snRNA-seq
%     specificData.transcriptomicData = readtable([dataFile filesep 'SCTnorm_control.csv']);
%       specificData.transcriptomicData = readtable([dataFile filesep 'scaled_control.csv']);
%     specificData.transcriptomicData = readtable([dataFile filesep 'DNcontrol_norexp1.csv']);
    % specificData.transcriptomicData = readtable([dataFile filesep 'DNcontrol_NB.csv']);
    specificData.transcriptomicData = readtable([dataFile filesep 'DNcontrol_NBnew.csv']);
    param.transcriptomicThreshold=-3;
    
end
specificData.transcriptomicData.expVal=specificData.transcriptomicData.expVal;
%% use proteomic data from Sarah Plum et al. 2020 (do not use because it is not suitbale to differentiate SYN and ASYN)
% ProteinFile='~/drive/bioenergeticsPD/fromXi/data/Proteomics/PMID: 33276480';
% 
% if any(strfind(bibliomicData,'PD')) & ~any(strfind(bibliomicData,'astro'))
%     specificData.proteomicData=readtable([ProteinFile filesep 'PDprotein.csv']);
%     param.thresholdP=22;%22 20 
%     
% elseif  ~any(strfind(bibliomicData,'PD')) & ~any(strfind(bibliomicData,'astro'))
%     specificData.proteomicData=readtable([ProteinFile filesep 'controlProtein.csv']);
%     param.thresholdP=20;%21 18
% end

%% 3.Technical parameters
%choose solver
solver='gurobi';
solverOK = changeCobraSolver(solver,'LP');
% parameters
param.TolMinBoundary = -1e3;
param.TolMaxBoundary =  1e3;
feasTol = getCobraSolverParams('LP', 'feasTol');
param.boundPrecisionLimit = feasTol * 100; % fluxEpsilon=1e-4;
param.closeIons = false; %not specified
% param.closeUptakes = false; %not specified
param.closeUptakes = true; % open the exchange rxns that same as the iDopa input data
% param.nonCoreSinksDemands = 'closeNone';
param.nonCoreSinksDemands = 'closeAll';
param.activeGenesApproach = 'oneRxnPerActiveGene';
% param.tissueSpecificSolver = 'fastCore'; % Choose fastcore
param.tissueSpecificSolver = 'thermoKernel'; % Choose thermoKernel
param.modelExtractionAlgorithm = 'thermoKernel'; % Choose thermoKernel
param.fluxEpsilon = feasTol * 100;% fluxEpsilon=1e-4;
param.fluxCCmethod = 'fastcc';
param.addCoupledRxns = 1; %  have been added in the preconditional recon model
param.thermoFluxEpsilon=param.fluxEpsilon;
param.weightsFromOmics=true;
param.curationOverOmics=true;
param.activeOverInactive=true;
param.inactiveGenesTranscriptomics = 1; 

%% 4.XomicsToModel pipeline to generate models
modelname=regexprep(bibliomicData,['_(\w+)ta'],'');
modelname=regexprep(modelname,'.xlsx','');
modelname=[modelname '_model'];

switch param.tissueSpecificSolver
    case 'fastCore'
resultsFolder = ['~/drive/bioenergeticsPD/fromXi/results/model/fastcore' filesep modelname];
    case 'thermoKernel'
resultsFolder = ['~/drive/bioenergeticsPD/fromXi/results/model/thermokernel' filesep modelname];
end

if ~isfolder(resultsFolder) 
    mkdir(resultsFolder);
end
cd(resultsFolder);

param.printLevel = 1;
param.debug = true;
if isunix()
    name = getenv('USER');
else
    name = getenv('username');
end
param.diaryFilename = [resultsFolder filesep datestr(now,30) '_' name '_diary.txt'];
%%
% Start the timer
startime=datetime('now');
[GeneratedModel, modelGenerationReport] = XomicsToModel(model, specificData, param);
% GeneratedModel.lb(contains(GeneratedModel.rxns,'DM_clpn_hs[c]'))=0.000649;
save([resultsFolder filesep [modelname '.mat']],'GeneratedModel')
% Stop the timer
endtime=datetime('now');
elapsedTime = endtime-startime;
