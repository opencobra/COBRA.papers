% This file describes step by step how to run the DEMETER pipeline on a
% single reconstruction. Adapt to your own organism as necessary.
% This script is intended for refinement of single reconstructions only.
% If you want to refine multiple draft reconstructions in batch, use the
% script "demeter_run_reconstructions_batch".
%
% Almut Heinken, 03/21

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

global CBTDIR

%% Step 1: preparing input data
% define the path to the draft reconstruction
draftFolder = [CBTDIR filesep 'papers' filesep '2021_demeter' filesep 'exampleDraftReconstructions'];
modelName='Lactobacillus_brevis_EW.RAST.fbamodel.sbml';

% define path to file with taxonomic information (create as needed)
infoFilePath = [CBTDIR filesep 'papers' filesep '2021_demeter' filesep 'example_infoFile.xlsx'];
infoFile = readtable(infoFilePath, 'ReadVariableNames', false);
infoFile = table2cell(infoFile);

% propagate experimental data
[infoFilePath,inputDataFolder] = prepareInputData(infoFilePath);

%% Step 2: refinement pipeline
% create an appropriate ID for the model
microbeID=adaptDraftModelID(modelName);

% load the model
draftModel = readCbModel([draftFolder filesep modelName]);

% create the model
[model,summary]=refinementPipeline(draftModel,microbeID, infoFilePath, inputDataFolder);

% save the model
save([microbeID '_refined'],'model');

% export the model as SBML file
writeCbModel(model, 'format', 'sbml', 'fileName', [microbeID '_revised']);

% save the summary of performed refinement
save([microbeID '_summary'],'summary');

% save the translated draft model for comparison
if contains(modelName,'sbml') || contains(modelName,'xml')
    draftModel = translateDraftReconstruction(draftModel);
    save([microbeID '_translatedDraft'],'draftModel');
end

%% Step 3: testing
%% Test of growth rates
% test aerobic and anaerobic growth on unlimited medium (consisting of 
% every compound the model can transport) and on a complex medium. The
% model should be able to grow anaerobically on the complex medium.
biomassReaction = model.rxns{find(strncmp(model.rxns,'bio',3)),1};
[AerobicGrowth, AnaerobicGrowth] = testGrowth(model, biomassReaction);

%% Test of ATP production
% Test ATP production under aerobic and anaerobic conditions on a complex
% medium. ATP production should not be higher than 150 mmol/g dry weight/hr
% aerobically and 100 mmol/g dry weight/hr.
[atpFluxAerobic, atpFluxAnaerobic] = testATP(model);

%% Testing of model predictions against experimental and comparativs data
% run the test suite on the model
testResultsRefined = runTestsOnModel(model, microbeID, inputDataFolder);
% remove all fields with no cases
fields = fieldnames(testResultsRefined);
lengths = structfun(@(x) size(x,2), testResultsRefined);
testResultsRefined = rmfield(testResultsRefined,fields(lengths<2));
% You will see that there are no false negative predictions in the test
% example after refinement.

% print a report on the test results as a PDF file. 
ncbiCol=find(strcmp(infoFile(1,:),'NCBI Taxonomy ID'));
if ~isempty(ncbiCol)
    ncbiID = infoFile(find(strcmp(infoFile(:,1),microbeID)),ncbiCol);
else
    ncbiID='';
end
outputFile = reportPDF(model, microbeID, biomassReaction, inputDataFolder, pwd, ncbiID);

% test the draft model for comparison. You will see that for the test
% example, there are multiple cases of false negative predictions.
testResultsDraft = runTestsOnModel(draftModel, microbeID, inputDataFolder);
% remove all fields with no cases
fields = fieldnames(testResultsDraft);
lengths = structfun(@(x) size(x,2), testResultsDraft);
testResultsDraft = rmfield(testResultsDraft,fields(lengths<2));

%% Step 4: debugging
% run additional gapfilling and elimination of futile cycles if neccessary
fields=fieldnames(testResultsRefined);

if atpFluxAerobic > 150 || atpFluxAnaerobic > 100 || AnaerobicGrowth(1,2) < 0.000001 || any(contains(fields,'FalseNegatives'))
    [revisedModel,gapfilledReactions,replacedReactions]=debugModel(model,testResultsRefined,inputDataFolder,microbeID,biomassReaction);
    
    % repeat tests on the revised model
    testResultsRevised = runTestsOnModel(revisedModel, microbeID, inputDataFolder);
    
    % export the revised model
    model=revisedModel;
    save([microbeID '_revised'],'model');
end

%% Step 5: computation of model properties
% Here, some features of the resulting reconstruction will be computed.
%% compute the uptake and secretion potential of the microbe
exRxns=model.rxns(find(strncmp(model.rxns,'EX_',3)));
% open all exchanges
model = changeRxnBounds(model, exRxns, -1000, 'l');
model = changeRxnBounds(model, exRxns, 1000, 'u');
try
    [minFlux, maxFlux, ~, ~] = fastFVA(model, 0, 'max', 'ibm_cplex', ...
        exRxns, 'S');
catch
    warning('fastFVA could not run, so fluxVariability is instead used. Consider installing fastFVA for shorter computation times.');
    cd(currentDir)
    [minFlux, maxFlux] = fluxVariability(model, 0, 'max', exRxns);
end

UptakeSecretion(:,1) = exRxns;
UptakeSecretion(:,2) = cellstr(num2str(minFlux));
UptakeSecretion(:,3) = cellstr(num2str(maxFlux));
% The variable "UptakeSecretion" contains the computed fluxes in 
% mmol/g dry weight/hr through exchange reactions in the model.

%% compute the internal metabolite biosynthesis potential
% get all internal metabolites
Metabolites=strrep(model.mets,'[c]','');
Metabolites=strrep(Metabolites,'[e]','');
Metabolites=unique(Metabolites);

modelTest=model;
% add sink reactions for each metabolite and optimize it
for k=1:length(Metabolites)
    if ~isempty(find(strcmp(modelTest.mets,[Metabolites{k} '[c]'])))
        if ~contains(modelTest.rxns,['DM_' Metabolites{k} '[c]'])
            modelTest=addDemandReaction(modelTest,[Metabolites{k} '[c]']);
        end
        modelTest=changeObjective(modelTest,['DM_' Metabolites{k} '[c]']);
        % do not count uptake of the metabolite
        modelTest=changeRxnBounds(modelTest,['EX_' Metabolites{k} '(e)'],0,'l');
        
        FBA=optimizeCbModel(modelTest,'max');
        Metabolites{k,2}=FBA.f;
    else
        Metabolites{k,2}=0;
    end
end
% The variable "Metabolites" contains the internal production potential in
% mmol/g dry weight/hr for each cytosolic metabolite present in the model.

