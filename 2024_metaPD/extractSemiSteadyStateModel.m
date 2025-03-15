function [extractModel, extractModelMetBool, extractModelRxnBool] = extractSemiSteadyStateModel(model,rxnWeights, metWeights, param)
%
% USAGE:
%   [thermoModel, thermoModelMetBool, thermoModelRxnBool] = extractSemiSteadyStateModel(model,rxnWeights, metWeights, param)
%
% INPUTS:
%  model.S:         m x n stoichiometric matrix       
%  model.ub:              
%  model.rxns:          
%  model.mets:          
%  model.lb:              
%  rxnWeights:          
%  metWeights:   
%
% OPTIONAL INPUTS
%  param.plotThermoKernelStats:
%  param.plotThermoKernelWeights:
%
% OUTPUTS:
%  extractModel:      
%  extractModelMetBool: m x 1 boolean, true where metabolite is part of extracted model 
%  extractModelRxnBool: n x 1 boolean, true where reaction is part of extracted model
%
% EXAMPLE:
%
% NOTE:
%
% Author(s):

if ~exist('param','var')
    param = struct();
end
if ~isfield(param,'plotThermoKernelStats')
    param.plotThermoKernelStats=1;
end
if ~isfield(param,'plotThermoKernelWeights')
    param.plotThermoKernelWeights=1;
end

[nMet,nRxns]=size(model.S);
if length(rxnWeights)~=nRxns
    error('rxnWeights must be same dimension as the columns of model.S')
end
if length(metWeights)~=nMet
    error('metWeights must be same dimension as the rows of model.S')
end
modelOri = model;

% identify the set of metabolites to relax steady state constraints by addition of reversible exchange reactions
relaxedMetBool = isnan(metWeights);

defaultRxnBound = max(abs(model.ub));
for i=1:nMet
    if relaxedMetBool(i)
        bool = strcmp(model.rxns,['EX_' model.mets{i}]) | strcmp(model.rxns,['EX_' model.mets{i} '[e]']);
        if any(bool)
            %open up the exchange
            model.lb(bool)=-defaultRxnBound;
            model.ub(bool)= defaultRxnBound;
            %no need for additional exchange reaction
            relaxedMetBool(i)=0;
        end
        
        bool = strcmp(model.rxns,['DM_' model.mets{i}]) | strcmp(model.rxns,['DM_' model.mets{i} '[e]']);
        if any(bool)
            %open up the demand
            model.ub(bool)= defaultRxnBound;
            %no need for additional exchange reaction
            relaxedMetBool(i)=0;
        end
        
        bool = strcmp(model.rxns,['sink_' model.mets{i}]) | strcmp(model.rxns,['sink_' model.mets{i} '[e]']);
        if any(bool)
            %open up the sink
            model.ub(bool)= defaultRxnBound;
            %no need for additional exchange reaction
            relaxedMetBool(i)=0;
        end
        
    end
end

    
lb = -defaultRxnBound*ones(nnz(relaxedMetBool),1);
ub =  defaultRxnBound*ones(nnz(relaxedMetBool),1);
[model, AddedExchRxn] = addExchangeRxn(model, model.mets(relaxedMetBool), lb, ub);
rxnWeights = [rxnWeights; zeros(length(AddedExchRxn),1)];

metWeights(isnan(metWeights))=0;

activeInactiveRxn=[];
presentAbsentMet=[];
[extractModel, extractModelMetBool, extractModelRxnBool] = ...
    thermoKernel(model, activeInactiveRxn, rxnWeights, presentAbsentMet, metWeights, param);