%% Prediction of perturbations to dopaminergic neuronal metabolism
% Load the model

if ~exist('modelName','var')
    modelName='iDopaNeuroC';
end

[~, ~] = changeCobraSolver('mosek', 'all', 0);
tcbmParam.solver = 'mosek';

if contains(char(java.lang.System.getProperty('user.name')),'aga')
    dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'exoMetDN' ...
        filesep 'data' filesep 'xomics'];
    modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'exoMetDN' ...
        filesep 'results' filesep 'codeResults' filesep 'iDN1'];
    resultsDir=['~' filesep 'work' filesep 'sbgCloud' filesep 'exoMetDN' ...
        filesep 'results' filesep 'codeResults' filesep 'iDN1' filesep 'iDopaNeuroCvsCT'];
else
    dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep ...
        'programReconstruction' filesep 'projects' filesep 'exoMetDN' ...
        filesep 'data' filesep 'xomics'];
    modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep ...
        'programReconstruction' filesep 'projects' filesep 'exoMetDN' ...
        filesep 'results' filesep 'codeResults' filesep 'iDN1'];
    resultsDir=['~' filesep 'work' filesep 'sbgCloud' filesep ...
        'programReconstruction' filesep 'projects' filesep 'exoMetDN' ...
        filesep 'results' filesep 'codeResults' filesep 'iDN1'];
end

%filenames
bibliomicData = 'bibliomicData.xlsx';
exometabolomicData = 'exometabolomicData.txt';
transcriptomicData = 'transcriptomicData.txt';
switch modelName
    case 'iDopaNeuroC'
        load([modelsDir filesep 'iDopaNeuroC' filesep 'iDopaNeuroC.mat']);
        resultsDir=fileparts([modelsDir filesep 'iDopaNeuroC' filesep 'iDopaNeuroC.mat']);
        model = iDopaNeuroC;
    case 'iDopaNeuroCT'
        load([modelsDir filesep 'iDopaNeuroCT' filesep 'iDopaNeuroCT.mat']);
        resultsDir=fileparts([modelsDir filesep 'iDopaNeuroCT' filesep 'iDopaNeuroT.mat']);
        model = iDopaNeuroCT;
end

if isfield(model, 'g0')
    model = rmfield(model, 'g0');
end
if isfield(model, 'g1')
    model = rmfield(model, 'g1');
end
%% 
% Save the vanilla model prior to any adjustments

model = changeRxnBounds(model, 'EX_gal[e]', 0, 'b'); %galactose not present in the medium
%modelOrig = model;
modelOrg = model;
%% 
% Manual correction of some exchange reactions lower boundaries

rxnsList = {'EX_atp[e]', 'EX_dopa[e]', 'EX_dopasf[e]', 'EX_h2o2[e]', 'EX_dopa4sf[e]', ...
    'EX_dopa4glcur[e]', 'EX_dopa3glcur[e]', 'EX_4glu56dihdind[e]', 'EX_5cysdopa[e]', 'EX_CE5025[e]'};
model = changeRxnBounds(model, rxnsList, 0, 'l');
modelOrig = model;
%% Read in metabolomic data

exoMet = readtable([dataFolder filesep 'exometabolomicData.txt']);
if isvar(exoMet,'rxnID')
    exoMet.Properties.VariableNames{'rxnID'} = 'rxns';
end
if isvar(exoMet,'metID')
    exoMet.Properties.VariableNames{'metID'} = 'mets';
end
if isvar(exoMet,'metName')
    exoMet.Properties.VariableNames{'metName'} = 'metNames';
end
%% 
% Hack to replace 0 0 with NaN

[mlt,nlt]=size(exoMet);
for i=1:mlt
    for j = 7:2:nlt
        if exoMet{i,j}==0 && exoMet{i,j+1}==0
            exoMet{i,j}=NaN;
            exoMet{i,j+1}=NaN;
        end
    end
end
%% Map the perturbation data onto the model

%[Bout,LIBkey,LOCAkey] = mapAontoB(Akey,Bkey,Ain,Bin)
[exoMet,LIBkey,LOCAkey] = mapAontoB(exoMet.rxns,model.rxns,exoMet);
exoMet = exoMet(LIBkey,:);
%% Extract the exometabolomic data
% Control exometabolomic data with glucose uptake

glcValidationData = exoMet(:, 1:4);
glcValidationData(:, end+1) = table(exoMet.mean);
glcValidationData(:, end+1) = table(exoMet.SD);
glcValidationData.Properties.VariableNames{'Var5'} = 'mean';
glcValidationData.Properties.VariableNames{'Var6'} = 'SD';
glcValidationData.Properties.Description='Exometabolomic Glucose medium';
glcValidationData = glcValidationData(~isnan(glcValidationData.mean),:);
glcValidationData = sortrows(glcValidationData, 'rxns');
%% 
% Galactose instead of glucose uptake

galValidationData = exoMet(:, 1:4);
galValidationData(:, end+1) = table(exoMet.mean_galactose);
galValidationData(:, end+1) = table(exoMet.sd_galactose);
galValidationData.Properties.VariableNames{'Var5'} = 'mean';
galValidationData.Properties.VariableNames{'Var6'} = 'SD';
galValidationData.Properties.Description='Exometabolomic Galactose medium';
galValidationData = galValidationData(~isnan(galValidationData.mean),:);
galValidationData = sortrows(galValidationData, 'rxns');
%% 
% Rotenone inhibition of Mitochondrial complex I

C1ValidationData = exoMet(:, 1:4);
C1ValidationData(:, end+1) = table(exoMet.mean_rotenone);
C1ValidationData(:, end+1) = table(exoMet.sd_rotenone);
C1ValidationData.Properties.VariableNames{'Var5'} = 'mean';
C1ValidationData.Properties.VariableNames{'Var6'} = 'SD';
C1ValidationData.Properties.Description='Exometabolomic C1 inhibition';
C1ValidationData = C1ValidationData(~isnan(C1ValidationData.mean),:);
C1ValidationData = sortrows(C1ValidationData, 'rxns');
%% 
% Oligomycin inhibition of Mitochonrial complex 5

C5ValidationData = exoMet(:, 1:4);
C5ValidationData(:, end+1) = table(exoMet.mean_oligomycin);
C5ValidationData(:, end+1) = table(exoMet.sd_oligomycin);
C5ValidationData.Properties.VariableNames{'Var5'} = 'mean';
C5ValidationData.Properties.VariableNames{'Var6'} = 'SD';
C5ValidationData.Properties.Description='Exometabolomic C5 inhibition';
C5ValidationData = C5ValidationData(~isnan(C5ValidationData.mean),:);
C5ValidationData = sortrows(C5ValidationData, 'rxns');
%% Control predictions
% Predict fluxes with EFBA

tcbmParam.method = 'fluxes';
tcbmParam.printLevel = 0;
model.osenseStr = 'min';
model.cf = 0;
model.cr = 0;
model.g = 2;
model.u0 = 0;
model.f = 1;
       
% control flux
[sol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
glcEFBAPrediction = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
glcEFBAPrediction.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
glcEFBAPrediction.v = sol.v;
glcEFBAPrediction.Properties.Description='EFBA Prediction, Glucose medium';
%% 
% Analyse the vanilla weights on the dual variables to net flux box constraints 
% to see which exchange reaction bounds are active

figure
feasTol = getCobraSolverParams('LP','feasTol');
min(sol.z_v)
max(sol.z_v)
histogram(sol.z_v(sol.z_v~=0 & abs(sol.z_v)>feasTol & ~model.SConsistentRxnBool))
%%
substantialDualToBoxConstraintsBool= abs(sol.z_v)>1e-4 & ~model.SConsistentRxnBool & model.lb~=0 & model.ub~=0 & ismember(model.rxns,exoMet.rxns);
fluxData = [model.lb,sol.v,model.ub,sol.y_v];
%printFluxVector(model, fluxData, nonZeroFlag, excFlag, sortCol, fileName, headerRow, formulaFlag, gprFlag)
printFluxVector(model, fluxData, substantialDualToBoxConstraintsBool, 0, 0, [], [], 1, 0)
%printConstraints(model,-inf, inf, substantialDualToBoxConstraintsBool)
%% Plot comparison of prediction with experimental data

param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.modelBounds=1;
plotExperimentalvsPredictedExchange(glcValidationData, glcEFBAPrediction, [], [], param)
%% 
% Generate a control uptake and secretion constrained model

model = modelOrig;
model = changeRxnBounds(model, 'EX_gal[e]', 0, 'l'); %galactose present
[modelUpt,modelSec] = generateModelUptSec(model,glcValidationData,param);
%% 
% Uptake constrained model

[sol, ~] = entropicFluxBalanceAnalysis(modelUpt, tcbmParam);
modelUpt_QEFBAGlc = table(model.rxns, model.rxnFormulas, modelUpt.lb, modelUpt.ub);
modelUpt_QEFBAGlc.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
modelUpt_QEFBAGlc.v = sol.v;
modelUpt_QEFBAGlc.Properties.Description='QEFBA Prediction, Glucose medium, modelUpt';
%% 
% Secretion constrained model

modelSec = changeRxnBounds(modelSec, 'EX_glc_D[e]', 0, 'b'); %glucose not present
modelSec = changeRxnBounds(modelSec, 'EX_gal[e]', -10000, 'l'); %galactose present
[sol, ~] = entropicFluxBalanceAnalysis(modelSec, tcbmParam);
modelSec_QEFBAGlc = table(model.rxns, model.rxnFormulas, modelSec.lb, modelSec.ub);
modelSec_QEFBAGlc.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
modelSec_QEFBAGlc.v = sol.v;
modelSec_QEFBAGlc.Properties.Description='QEFBA Prediction, Glucose medium, modelSec';
%%
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.measuredBounds=1;
param.predictedBounds=0;
param.otherBounds=0;
param.extraBounds=0;
plotExperimentalvsPredictedExchange(glcValidationData, modelUpt_QEFBAGlc, modelSec_QEFBAGlc, [], param)
%% Fit the model exchanges to the exometabolomic constraints from the glucose medium condition
% Fit a control model to experimentally measured fluxes while maximising flux 
% entropy and maximising a linear objective on reaction flux. Add the additional 
% quadratic terms to the model to be solved with EFBA and relax the bounds on 
% the reactions with exometabolomic constraints from the control condition

param.alpha=1e5;
param.printLevel=1;
param.relaxBounds=1;
model = addExoMetToEFBA(modelOrig,glcValidationData,param);
%%
[sol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
glcQEFBAPrediction = glcEFBAPrediction;
glcQEFBAPrediction.v = sol.v;
glcQEFBAPrediction.Properties.Description='QEFBA Prediction, Glucose medium';
disp(sol)
posZw = sol.z_v>feasTol;
negZw = sol.z_v<feasTol;
disp(sum(sol.z_v(negZw)'*model.lb(negZw) + sol.z_v(posZw)'*model.ub(posZw)))
%%
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.modelBounds=1;
plotExperimentalvsPredictedExchange(glcValidationData, glcEFBAPrediction, glcQEFBAPrediction, [], param)
%% 
% Leave one out cross validation on control glucose uptake model

model = changeRxnBounds(modelOrig, 'EX_gal[e]', 0, 'b'); %galactose not present in the medium
LOOCVparam.alpha=1e5;
LOOCVparam.printLevel=1;
LOOCVparam.relaxBounds=1;
[glcLOOCV,~] = exoMetLOOCV(model,glcValidationData,LOOCVparam);
%%
glcQEFBAPrediction = glcEFBAPrediction;
glcQEFBAPrediction.v = glcLOOCV.QEFBA;
glcQEFBAPrediction.Properties.Description='QEFBA Prediction, Glucose medium';
glcLOOCVPrediction = glcEFBAPrediction;
glcLOOCVPrediction.v = glcLOOCV.LOOCV;
glcLOOCVPrediction.Properties.Description='LOOCV, QEFBA Prediction, Glucose medium';
plotExperimentalvsPredictedExchange(glcValidationData, glcEFBAPrediction, glcQEFBAPrediction, glcLOOCVPrediction, param)
%% Plot ATP balance in the control model
% Uses the param.treshold_v to display only the reactions that carry significant 
% flux

param.treshold_v = 1;
param.NodeLabels = 'rxnNames';
[graphGlc,summaryGraphGlc] = metUtilisation(model, 'atp[c]', glcQEFBAPrediction.v, 1, param);
param.treshold_v = 0.1;
[graphGlcNADH,summaryGraphGlcNADH] = metUtilisation(model, 'nadh[m]', glcQEFBAPrediction.v, 1, param);
%% Switch to uptake of galactose instead of glucose
% Check that quadratic penalisation to measured exchanges in galactose medium 
% is possible

param.alpha=1e5;
param.printLevel=0;
param.relaxBounds=1;
model = addExoMetToEFBA(modelOrig,galValidationData,param);
model = changeRxnBounds(model, 'EX_glc_D[e]', 0, 'b'); %glucose not present
model = changeRxnBounds(model, 'EX_gal[e]', -10000, 'l'); %galactose present
model.c(ismember(model.rxns,'ATPM'))=0;%1000;
[sol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
galQEFBAPrediction = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
galQEFBAPrediction.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
galQEFBAPrediction.v = sol.v;
galQEFBAPrediction.Properties.Description='QEFBA Prediction, Galactose medium';
disp(sol.stat)
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.measuredBounds=1;
param.predictedBounds=0;
param.otherBounds=1;
param.extraBounds=0;
plotExperimentalvsPredictedExchange(galValidationData, galQEFBAPrediction, [], [], param)
%% 
% Generate an uptake and secretion constrained model

model = modelOrig;
param.quadPenalty = 1;
param.alpha = 1e5;
[modelUpt,modelSec] = generateModelUptSec(model,galValidationData,param);
%% 
% Uptake constrained model

modelUpt = changeRxnBounds(modelUpt, 'EX_glc_D[e]', 0, 'b'); %glucose not present
modelUpt = changeRxnBounds(modelUpt, 'EX_gal[e]', -10000, 'l'); %galactose present
[sol, ~] = entropicFluxBalanceAnalysis(modelUpt, tcbmParam);
modelUpt_qeFBAGal = table(model.rxns, model.rxnFormulas, modelUpt.lb, modelUpt.ub);
modelUpt_qeFBAGal.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
modelUpt_qeFBAGal.v = sol.v;
modelUpt_qeFBAGal.Properties.Description='QEFBA Prediction, Galactose medium, modelUpt';
%% 
% Secretion constrained model

modelSec = changeRxnBounds(modelSec, 'EX_glc_D[e]', 0, 'b'); %glucose not present
modelSec = changeRxnBounds(modelSec, 'EX_gal[e]', -10000, 'l'); %galactose present
[sol, ~] = entropicFluxBalanceAnalysis(modelSec, tcbmParam);
modelSec_qeFBAGal = table(model.rxns, model.rxnFormulas, modelSec.lb, modelSec.ub);
modelSec_qeFBAGal.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
modelSec_qeFBAGal.v = sol.v;
modelSec_qeFBAGal.Properties.Description='QEFBA Prediction, Galactose medium, modelSec';
%%
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.measuredBounds=1;
param.predictedBounds=0;
param.otherBounds=0;
param.extraBounds=0;
plotExperimentalvsPredictedExchange(galValidationData, modelUpt_qeFBAGal, modelSec_qeFBAGal, [], param)
%% 
% Leave one out cross validation on control galactose uptake model

model = changeRxnBounds(modelOrig, 'EX_glc_D[e]', 0, 'b'); %glucose not present
model = changeRxnBounds(model, 'EX_gal[e]', -10000, 'l'); %galactose present
LOOCVparam.alpha=1e5;
LOOCVparam.printLevel=1;
LOOCVparam.relaxBounds=1;
[glcLOOCV,~] = exoMetLOOCV(model,galValidationData,LOOCVparam);
%%
galLOOCVPrediction = glcEFBAPrediction;
galLOOCVPrediction.v = glcLOOCV.LOOCV;
galLOOCVPrediction.Properties.Description='LOOCV, QEFBA Prediction, Galactose medium';
plotExperimentalvsPredictedExchange(galValidationData, galQEFBAPrediction, galLOOCVPrediction, [], param)
%% Plot ATP balance in the galactose model and compare to glucose model
% Uses the param.treshold_v to display only the reactions that carry significant 
% flux

param.treshold_v = 1;
[graphGal,summaryGraphGal] = metUtilisation(model, 'atp[c]', galQEFBAPrediction.v, 1, param);
[graphGluGal,summaryGraphGluGal] = metUtilisation(model, 'atp[c]', [glcQEFBAPrediction.v galQEFBAPrediction.v], 1, param);
%% Inhibition of mitochondrial complex I

if contains('NADH2_u10minew',modelOrig.rxns)
    NADH2_u10m = 'NADH2_u10minew'; % '5 h[m] + nadh[m] + q10[m]  -> 4 h[c] + nad[m] + q10h2[m] '
else
    NADH2_u10m = 'NADH2_u10m';
end
%% 
% Check that quadratic penalisation to measured exchanges with C1 inhibition 
% is possible

param.alpha=1e5;
param.printLevel=0;
model = changeRxnBounds(modelOrig, NADH2_u10m, 0, 'b'); %inhibit complex 1
model = changeRxnBounds(model, 'EX_gal[e]', 0, 'b'); %galactose not present
model = addExoMetToEFBA(model,C1ValidationData,param);
[sol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
C1QEFBAPrediction = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
C1QEFBAPrediction.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
C1QEFBAPrediction.v = sol.v;
C1QEFBAPrediction.Properties.Description='QEFBA Prediction, Complex 1 inhibition';
%disp(sol)
%%
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.measuredBounds=1;
param.predictedBounds=0;
param.otherBounds=1;
param.extraBounds=0;
plotExperimentalvsPredictedExchange(C1ValidationData, C1QEFBAPrediction, [], [], param)
%set(gcf,'Visible','on')
%% 
% Generate an uptake and secretion constrained model

param.quadPenalty = 1;
param.alpha = 1e5;
[modelUpt,modelSec] = generateModelUptSec(model,C1ValidationData,param);
%% 
% Uptake constrained model

[sol, ~] = entropicFluxBalanceAnalysis(modelUpt, tcbmParam);
modelUpt_QEFBAC1 = table(model.rxns, model.rxnFormulas, modelUpt.lb, modelUpt.ub);
modelUpt_QEFBAC1.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
modelUpt_QEFBAC1.v = sol.v;
modelUpt_QEFBAC1.Properties.Description='QEFBA Prediction, C1 inhibition, modelUpt';
%% 
% Secretion constrained model

[sol, ~] = entropicFluxBalanceAnalysis(modelSec, tcbmParam);
modelSec_QEFBAC1 = table(model.rxns, model.rxnFormulas, modelSec.lb, modelSec.ub);
modelSec_QEFBAC1.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
modelSec_QEFBAC1.v = sol.v;
modelSec_QEFBAC1.Properties.Description='QEFBA Prediction, C1 inhibition, modelSec';
%%
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.measuredBounds=1;
param.predictedBounds=0;
param.otherBounds=0;
param.extraBounds=0;
plotExperimentalvsPredictedExchange(C1ValidationData, modelUpt_QEFBAC1, modelSec_QEFBAC1, [], param)
%set(gcf,'Visible','on')
%% 
% Leave one out cross validation on Complex I inhibition model

model = changeRxnBounds(modelOrig, NADH2_u10m, 0, 'b'); %inhibit complex 1
model = changeRxnBounds(model, 'EX_gal[e]', 0, 'b'); %galactose not present
LOOCVparam.alpha=1e5;
LOOCVparam.printLevel=1;
LOOCVparam.relaxBounds=1;
[C1LOOCV,~] = exoMetLOOCV(model,C1ValidationData,LOOCVparam);
%%
C1LOOCVPrediction = glcEFBAPrediction;
C1LOOCVPrediction.v = glcLOOCV.LOOCV;
C1LOOCVPrediction.Properties.Description='LOOCV, QEFBA Prediction, C1 inhibiton';
plotExperimentalvsPredictedExchange(C1ValidationData, C1QEFBAPrediction, C1LOOCVPrediction, [], param)
%% Plot ATP balance in the Complex I inhibition model and compare to it to the glucose model
% Uses the param.treshold_v to display only the reactions that carry significant 
% flux

param.treshold_v = 1;
[graphC1,summaryGraphC1] = metUtilisation(model, 'atp[c]', C1QEFBAPrediction.v, 1, param);
[graphGluC1,summaryGraphGluC1] = metUtilisation(model, 'atp[c]', [glcQEFBAPrediction.v C1QEFBAPrediction.v], 1, param);
%% Inhibition of mitochondrial complex V

if contains('ATPS4minew',modelOrig.rxns)
    ATPS4m = 'ATPS4minew'; % '4 h[c] + adp[m] + pi[m]  -> h2o[m] + 3 h[m] + atp[m] '
    NADH2_u10m = 'NADH2_u10minew'; % '5 h[m] + nadh[m] + q10[m]  -> 4 h[c] + nad[m] + q10h2[m] '
else
    ATPS4m = 'ATPS4m'; %adp[m] + 4.0 h[c] + pi[m] -> atp[m] + 3.0 h[m] + h2o[m] original, not found in moded
    NADH2_u10m = 'NADH2_u10m';
end
%% 
% Check that quadratic penalisation to measured exchanges with C5 inhibition 
% is possible

param.alpha = 1e5;
param.printLevel=0;
model = changeRxnBounds(modelOrig, ATPS4m, 0, 'b'); %inhibit complex V
model = changeRxnBounds(model, 'EX_gal[e]', 0, 'b'); %galactose not present
model = addExoMetToEFBA(model,C5ValidationData,param);
[sol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
C5QEFBAPrediction = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
C5QEFBAPrediction.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
C5QEFBAPrediction.v = sol.v;
C5QEFBAPrediction.Properties.Description='QEFBA Prediction, Complex 5 inhibition';
%disp(sol)
%%
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.measuredBounds=1;
param.predictedBounds=0;
param.otherBounds=1;
param.extraBounds=0;
plotExperimentalvsPredictedExchange(C5ValidationData, C5QEFBAPrediction, [], [], param)
%set(gcf,'Visible','on')
%% 
% Generate an uptake and secretion constrained model

param.quadPenalty = 1;
param.alpha = 1e5;
model = changeRxnBounds(modelOrig, ATPS4m, 0, 'b'); %inhibit complex V
model = changeRxnBounds(model, 'EX_gal[e]', 0, 'b'); %galactose not present
[modelUpt,modelSec] = generateModelUptSec(model,C5ValidationData,param);
%% 
% Uptake constrained model

[sol, ~] = entropicFluxBalanceAnalysis(modelUpt, tcbmParam);
modelUpt_QEFBAC5 = table(model.rxns, model.rxnFormulas, modelUpt.lb, modelUpt.ub);
modelUpt_QEFBAC5.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
modelUpt_QEFBAC5.v = sol.v;
modelUpt_QEFBAC5.Properties.Description='QEFBA Prediction, C5 inhibition, modelUpt';
%% 
% Secretion constrained model

[sol, ~] = entropicFluxBalanceAnalysis(modelSec, tcbmParam);
modelSec_QEFBAC5 = table(model.rxns, model.rxnFormulas, modelSec.lb, modelSec.ub);
modelSec_QEFBAC5.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
modelSec_QEFBAC5.v = sol.v;
modelSec_QEFBAC5.Properties.Description='QEFBA Prediction, C5 inhibition, modelSec';
%%
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.measuredBounds=1;
param.predictedBounds=0;
param.otherBounds=0;
param.extraBounds=0;
plotExperimentalvsPredictedExchange(C5ValidationData, modelUpt_QEFBAC5, modelSec_QEFBAC5, [], param)
%set(gcf,'Visible','on')
%% 
% Leave one out cross validation on Complex V inhibition model

model = changeRxnBounds(modelOrig, ATPS4m, 0, 'b'); %inhibit complex V
model = changeRxnBounds(model, 'EX_gal[e]', 0, 'b'); %galactose not present
LOOCVparam.alpha=1e5;
LOOCVparam.printLevel=1;
LOOCVparam.relaxBounds=1;
[C5LOOCV,LIBkey] = exoMetLOOCV(model,C5ValidationData,LOOCVparam);
%%
C5LOOCVPrediction = glcEFBAPrediction;
C5LOOCVPrediction.v = C5LOOCV.LOOCV;
C5LOOCVPrediction.Properties.Description='LOOCV, QEFBA Prediction, Complex 5 inhibiton';
plotExperimentalvsPredictedExchange(C5ValidationData, C5QEFBAPrediction, C5LOOCVPrediction, [], param)
%% Plot ATP balance in the Complex V inhibition model and compare to it to the glucose model
% Uses the param.treshold_v to display only the reactions that carry significant 
% flux

param.treshold_v = 1;
[graphC5,summaryGraphC5] = metUtilisation(model, 'atp[c]', C5QEFBAPrediction.v, 1, param);
[graphGluC5,summaryGraphGluC5] = metUtilisation(model, 'atp[c]', [glcQEFBAPrediction.v C5QEFBAPrediction.v], 1, param);
%% Inhibition of GBA1 gene

param.alpha = 1e5;
param.printLevel=0;
QEFBA=1;
if QEFBA 
    model = addExoMetToEFBA(modelOrig,glcValidationData,param);
else
    model = modelOrig;
end
%% 
% Control prediction

[sol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
glcQEFBAPrediction = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
glcQEFBAPrediction.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
glcQEFBAPrediction.v = sol.v;
glcQEFBAPrediction.Properties.Description='QEFBA Prediction, glucose medium';
glcQEFBAPrediction.labels = model.rxnNames;
%% 
% Control prediction 2

[sol2, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
glcQEFBAPrediction2 = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
glcQEFBAPrediction2.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
glcQEFBAPrediction2.v = sol2.v;
glcQEFBAPrediction2.Properties.Description='QEFBA Prediction 2, glucose medium';
glcQEFBAPrediction2.labels = model.rxnNames;
%% 
% Check reproducibility of control prediction

param.rxns=model.rxns(abs(glcQEFBAPrediction.v-glcQEFBAPrediction2.v)>1);
plotPredictedExchange(glcQEFBAPrediction, glcQEFBAPrediction2, [], param)
%% 
% Identify GBA1 related reactions

gbaId = model.genes(find(~cellfun(@isempty,strfind(model.genes,'2629')))); 
[~, ~, GBArxns, ~] = deleteModelGenes(model, gbaId, 0);
printConstraints(model,[],[],GBArxns)
%% 
% 
% 
% Optimise for GBA1 related reactions

param.alpha = 1e5;
param.printLevel=0;
if QEFBA 
    model = addExoMetToEFBA(modelOrig,glcValidationData,param);
else
    model = modelOrig;
end

optGBA1QEFBAPrediction = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
optGBA1QEFBAPrediction.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
if 0
    model.c(ismember(model.rxns,GBArxns))=60;
    optGBA1QEFBAPrediction.Properties.Description='QEFBA Prediction, GBA maximisation';
else
    model = changeRxnBounds(model, 'GBA',10,'l');
    optGBA1QEFBAPrediction.Properties.Description='QEFBA Prediction, GBA lb = 10';
end
[sol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
optGBA1QEFBAPrediction.v = sol.v;
optGBA1QEFBAPrediction.v(ismember(model.rxns,'GBA'))
optGBA1QEFBAPrediction.v(ismember(model.rxns,'GBAl'))
figure;histogram(abs(optGBA1QEFBAPrediction.v-glcQEFBAPrediction.v))
param.saveFigures=0;
param.expOrder = 1; %1  = order reactions by value of experimental flux
param.predictedBounds=0;
param.otherBounds=1;
param.extraBounds=0;
param.rxns=model.rxns(abs(optGBA1QEFBAPrediction.v-glcQEFBAPrediction.v)>1);
plotPredictedExchange(glcQEFBAPrediction, optGBA1QEFBAPrediction, [], param)
%% 
% Quadratic penalisation to measured control exchanges with GBA1 inhibition

param.alpha = 1e5;
param.printLevel=0;
if QEFBA
    model = addExoMetToEFBA(modelOrig,glcValidationData,param);
else
    model = modelOrig;
end
%% 
% GBA deletion 100%

if 1
    model = changeRxnBounds(model, GBArxns,0,'b');
end
[sol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
GBA1delQEFBAPrediction = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
GBA1delQEFBAPrediction.Properties.VariableNames = {'rxns', 'rxnFormulas', 'lb', 'ub'};
GBA1delQEFBAPrediction.v = sol.v;
GBA1delQEFBAPrediction.Properties.Description='QEFBA Prediction, GBA inhibition';
figure;histogram(abs(glcQEFBAPrediction.v-GBA1delQEFBAPrediction.v))
%%
param.saveFigures=0;
param.predictedBounds=1;
param.otherBounds=0;
param.extraBounds=0;
param.exchangesOnly=0;
%bool = abs(optGBA1QEFBAPrediction.v-GBA1delQEFBAPrediction.v)./sum([optGBA1QEFBAPrediction.v,GBA1delQEFBAPrediction.v],2)>0.5;
bool = abs(optGBA1QEFBAPrediction.v-GBA1delQEFBAPrediction.v)>3;
nnz(bool)
nnz(abs(optGBA1QEFBAPrediction.v-GBA1delQEFBAPrediction.v)>0.5) - nnz(bool)
param.rxns=[model.rxns(bool);'GBA'];
optGBA1QEFBAPrediction.labels = model.rxnNames;
plotPredictedExchange(optGBA1QEFBAPrediction, GBA1delQEFBAPrediction, [], param)
%% Plot ATP balance in the GBA inhibition model and compare to it to the activated GBA model
% Uses the param.treshold_v to display only the reactions that carry significant 
% flux

param.treshold_v = 1;
[graphGBA,summaryGraphGBA] = metUtilisation(model, 'atp[c]', GBA1delQEFBAPrediction.v, 1, param);
[graphGluGBA,summaryGraphGluGBA] = metUtilisation(model, 'atp[c]', [optGBA1QEFBAPrediction.v GBA1delQEFBAPrediction.v], 1, param);
%% Compare reaction fluxes in various perturbed conditions

%% Plot specific intracellular changes
figure()
% Top two plots
H = tiledlayout(2,2);
H.TileSpacing = 'compact';
H.Padding = 'tight';
nexttile

Xlabels = [model.rxnNames(ismember(model.rxns,'PGK')),model.rxnNames(ismember(model.rxns,'PYK'))];
V = zeros(4,2);
V(1,1) = glcQEFBAPrediction.v(ismember(model.rxns,'PGK'));
V(2,1) = C1QEFBAPrediction.v(ismember(model.rxns,'PGK'));
V(3,1) = C5QEFBAPrediction.v(ismember(model.rxns,'PGK'));
V(4,1) = galQEFBAPrediction.v(ismember(model.rxns,'PGK'));
V(1,2) = glcQEFBAPrediction.v(ismember(model.rxns,'PYK'));
V(2,2) = C1QEFBAPrediction.v(ismember(model.rxns,'PYK'));
V(3,2) = C5QEFBAPrediction.v(ismember(model.rxns,'PYK'));
V(4,2) = galQEFBAPrediction.v(ismember(model.rxns,'PYK'));

bar(1:2,abs(V), 1)
ylim([0 max(max(abs(V)))+10])
% Add title and axis labels
title('A. Glycolysis','FontSize',14)
set(gca, 'XTick', 1:2, 'XTickLabel', Xlabels,'FontSize',14)
ylabel('Flux (uMol/gDW/hr)','FontSize',14)
% Add a legend
legend('Glucose media', 'Complex I inhibition', 'Complex V inhibition', 'Galactose media', 'Location', 'northeast')

nexttile

Xlabels = [model.rxnNames(ismember(model.rxns,'PGL')),model.rxnNames(ismember(model.rxns,'GND'))];
V = zeros(4,2);
V(1,1) = glcQEFBAPrediction.v(ismember(model.rxns,'PGL'));
V(2,1) = C1QEFBAPrediction.v(ismember(model.rxns,'PGL'));
V(3,1) = C5QEFBAPrediction.v(ismember(model.rxns,'PGL'));
V(4,1) = galQEFBAPrediction.v(ismember(model.rxns,'PGL'));
V(1,2) = glcQEFBAPrediction.v(ismember(model.rxns,'GND'));
V(2,2) = C1QEFBAPrediction.v(ismember(model.rxns,'GND'));
V(3,2) = C5QEFBAPrediction.v(ismember(model.rxns,'GND'));
V(4,2) = galQEFBAPrediction.v(ismember(model.rxns,'GND'));

bar(1:2, V, 1)
ylim([0 max(max(abs(V)))+1])
%box off
set(gca, 'XTick', 1:2, 'XTickLabel', Xlabels,'FontSize',14)
% Add title and axis labels
%title('Predicted flux of the core energy-related reactions in control and complex I inhibition.')
ylabel('Flux (uMol/gDW/hr)','FontSize',14)
title('B. Pentose Phosphate Pathway','FontSize',14)


    
% Plot that spans
nexttile([1 2])
% ATPS4m = 'ATPS4minew'; % '4 h[c] + adp[m] + pi[m]  -> h2o[m] + 3 h[m] + atp[m] '
% FADH2ETC = Complex Ii Reaction for Respiratory Chain fadh2[m] + q10[m] -> fad[m] + q10h2[m] 
% NADH2_u10m = 'NADH2_u10minew'; % '5 h[m] + nadh[m] + q10[m]  -> 4 h[c] + nad[m] + q10h2[m] '
%Xlabels = [model.rxnNames(ismember(model.rxns,NADH2_u10m)),model.rxnNames(ismember(model.rxns,'FADH2ETC')),model.rxnNames(ismember(model.rxns,ATPS4m))];
Xlabels = {'Complex I','Complex II','Complex III','Complex IV','Complex V'};
V = zeros(4,3);

V(1,1) = glcQEFBAPrediction.v(ismember(model.rxns,NADH2_u10m));
V(2,1) = C1QEFBAPrediction.v(ismember(model.rxns,NADH2_u10m));
V(3,1) = C5QEFBAPrediction.v(ismember(model.rxns,NADH2_u10m));
V(4,1) = galQEFBAPrediction.v(ismember(model.rxns,NADH2_u10m));
V(1,2) = glcQEFBAPrediction.v(ismember(model.rxns,'FADH2ETC'));
V(2,2) = C1QEFBAPrediction.v(ismember(model.rxns,'FADH2ETC'));
V(3,2) = C5QEFBAPrediction.v(ismember(model.rxns,'FADH2ETC'));
V(4,2) = galQEFBAPrediction.v(ismember(model.rxns,'FADH2ETC'));
V(1,3) = glcQEFBAPrediction.v(ismember(model.rxns,'CYOR_u10minew'));
V(2,3) = C1QEFBAPrediction.v(ismember(model.rxns,'CYOR_u10minew'));
V(3,3) = C5QEFBAPrediction.v(ismember(model.rxns,'CYOR_u10minew'));
V(4,3) = galQEFBAPrediction.v(ismember(model.rxns,'CYOR_u10minew'));
V(1,4) = glcQEFBAPrediction.v(ismember(model.rxns,'CYOOm3inew'))+glcQEFBAPrediction.v(ismember(model.rxns,'CYOOm2inew'));
V(2,4) = C1QEFBAPrediction.v(ismember(model.rxns,'CYOOm3inew'))+C1QEFBAPrediction.v(ismember(model.rxns,'CYOOm2inew'));
V(3,4) = C5QEFBAPrediction.v(ismember(model.rxns,'CYOOm3inew'))+C5QEFBAPrediction.v(ismember(model.rxns,'CYOOm2inew'));
V(4,4) = galQEFBAPrediction.v(ismember(model.rxns,'CYOOm3inew'))+galQEFBAPrediction.v(ismember(model.rxns,'CYOOm2inew'));
V(1,5) = glcQEFBAPrediction.v(ismember(model.rxns,ATPS4m));
V(2,5) = C1QEFBAPrediction.v(ismember(model.rxns,ATPS4m));
V(3,5) = C5QEFBAPrediction.v(ismember(model.rxns,ATPS4m));
V(4,5) = galQEFBAPrediction.v(ismember(model.rxns,ATPS4m));

bar(1:5, abs(V), 1)

ylim([0 max(max(abs(V)))+1])
%box off
set(gca, 'XTick', 1:5, 'XTickLabel', Xlabels,'FontSize',14)
% Add title and axis labels
%title('Predicted flux of the core energy-related reactions in control and complex I inhibition.')
ylabel('Flux (uMol/gDW/hr)','FontSize',14)
% Add a legend
title('C. Oxidative phosphorylation','FontSize',14)
%iDopaNeuroC_coreRxns_perturbations.png
%% Summary table with predictions of fluxes for all perturbations

summaryTableAll = table(modelOrig.rxns, modelOrig.rxnNames, printRxnFormula(modelOrig, modelOrig.rxns), ...
    modelOrig.subSystems, modelOrig.lb, modelOrig.ub, glcQEFBAPrediction.v, galQEFBAPrediction.v, ...
    C1QEFBAPrediction.v, C5QEFBAPrediction.v, optGBA1QEFBAPrediction.v, GBA1delQEFBAPrediction.v);
summaryTableAll.Properties.VariableNames = {'rxns', 'rxnNames', 'rxnFormulas', 'subSystems', 'lb', 'ub', ...
    [modelName ' Control'], [modelName ' Galactose'], [modelName ' Complex I inhibition'], ...
    [modelName ' Complex V inhibition'], [modelName ' GBA1 activation'], [modelName ' GBA1 inhibition']};

if contains(char(java.lang.System.getProperty('user.name')),'aga')
    writetable(summaryTableAll,'~/work/sbgCloud/exoMetDN/papers/v22/SM/SM2.xlsx','Sheet',[modelName ' Flux Estimation']);
else
    writetable(summaryTableAll,[resultsDir filesep 'SM2.xlsx'],'Sheet',[modelName ' Flux Estimation']);
end