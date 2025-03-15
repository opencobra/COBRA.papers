function [CharacteristicTable] = ModelCharacteristic(multimodels,printFlag)
%To generate a table with the Characteristic results for each model
% USAGE:
%
%       [CharacteristicTable] = ModelCharacteristic(models,printFlag)
%
% INPUT:
%       multiModels:     struct format with models that need to compare
%                        e.g.multimodels.model1
%                           multimodels.model2 ...
%       printFlag:       1 if information should be printed to a table.
%                        Default = 0
%
% OUTPUT
%       CharacteristicTable: a table that contains the basic number of
%       rxns, mets and genes of each model
%
%EXAMPLE:
%
%    [CharacteristicTable] = odelCharacteristic(multimodels,1)
%
%NOTE:
%
%    This function is used to explore SYN/ASYN models from xomicsToModel
%    pipeline
%
%Author(s): - Xi Luo

if ~exist('printFlag', 'var') || isempty(printFlag)
    printFlag = 0;
elseif (~isnumeric(printFlag) & ~islogical(printFlag))
    error('printFlag should be a number or a bool')
end

Newmodels=multimodels;

if isstruct(Newmodels)
    % model characteristics
    Characteristic=cell(10,length(fieldnames(Newmodels)));
    modelname=fieldnames(Newmodels);
    for i=1:length(modelname)
        model=Newmodels.(modelname{i});
        % Total num of genes
        Characteristic{1,i}=length(unique(model.genes));
        % Total num of genes in mito
        Mrxns=findRxnFromCompartment(model,'[m]');
        Mgenes=findGenesFromRxns(model,Mrxns(:,1));
        Mgenes(find(cellfun(@isempty,Mgenes)))=[];
        uniquegenes = [];
        for m = 1:length(Mgenes)
            uniquegenes = [uniquegenes; Mgenes{m}(:)];
        end
        Characteristic{2,i}=length(unique(uniquegenes));
        % Total num of rxns
        Characteristic{3,i}=length(unique(model.rxns));
        % Total num of rxns in mito
        Characteristic{4,i}=size(Mrxns,1);
        % Total num of unique mets
        Characteristic{5,i}=length(unique(regexprep(model.mets,'(\[\w\])','')));
        % Total num of unique mets in mito
        compartments=regexp(model.mets, '\[(.*?)\]', 'match');
        Characteristic{6,i}=sum(string(compartments)=='[m]');
        % Total num of transport rxns
        Trans={'Transport, mitochondrial','Transport, extracellular','Transport, endoplasmic reticular'...
            'Transport, golgi apparatus','Transport, lysosomal','Transport, nuclear','Transport, peroxisomal'};
        TransRxns=[];
        for j=1:length(Trans)
            rxns=findRxnsFromSubSystem(model,Trans{j});
            TransRxns=[TransRxns ; rxns];
        end
        Characteristic{7,i}=length(unique(TransRxns));
        % Total num of exchange rxns
        if ~isfield(model,'ExchRxnBool') || ~isfield(model,'DMRxnBool') || ~isfield(model,'SinkRxnBool')
            model=findSExRxnInd(model);
        end
%         [selExc,selUpt] = findExcRxns(ComplexImodel);
%         indSyn=find(selUpt);
        Characteristic{8,i}=length(unique(model.rxns(model.ExchRxnBool==1)));
        % Characteristic{9,i}=length(unique(model.rxns(model.DMRxnBool==1)));
        %Characteristic{10,i}=length(unique(model.rxns(model.SinkRxnBool==1)));  
        % Total num of subsystems
        Characteristic{9,i}=length(unique(string(model.subSystems)));
        % Total num of compartments
        Characteristic{10,i}=length(unique(string(compartments)));
    end

Characteristic=cell2table(Characteristic);
Characteristic.Properties.VariableNames=fieldnames(Newmodels);
Characteristic.Properties.RowNames={'Total num of genes','Total num of genes in the mito','Total num of reactions',...
    'Total num of reactions in the mito','Total num of metabolites', 'Total num of metabolites in the mito',...
    'Total num of transport reactions','Total num of exchange reactions','Total num of subsystems','Total num of Compartments'};
CharacteristicTable=Characteristic;
else
    disp('please check the input model')
end
%% Print tables with output if printFlag = 1
if printFlag ==1
    CharacteristicTable
end
