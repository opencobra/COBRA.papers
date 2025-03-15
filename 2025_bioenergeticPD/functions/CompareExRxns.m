function [compareExRxnsTable] = CompareExRxns(Allmodels)
%UNTITLED3 Summary of this function goes here
if ~exist('Allmodels','var') || ~isstruct(Allmodels)
    error('check the input models')
end

if any(contains(fieldnames(Allmodels),{'SYNPD'}))
    names = [fieldnames(Allmodels.SYN); fieldnames(Allmodels.SYNPD)];
    Allmodels.allSYN = cell2struct([struct2cell(Allmodels.SYN); struct2cell(Allmodels.SYNPD)], names, 1);
    names = [fieldnames(Allmodels.ASYN); fieldnames(Allmodels.ASYNPD)];
    Allmodels.allASYN = cell2struct([struct2cell(Allmodels.ASYN); struct2cell(Allmodels.ASYNPD)], names, 1); 
    Allmodels=rmfield(Allmodels,{'SYN','SYNPD','ASYN','ASYNPD'});
end


type=fieldnames(Allmodels);
type=type(~contains(type,'Astro'));% ignore astro models
for i =1:length(type)
    models=fieldnames(Allmodels.(type{i}));
    models=models(~contains(models,'constrain'));% ignore unconstrained models
        % calculate for old and new models
        % SYN1 vs SYN2
        group.name1=models(~contains(models,'PD'));
        % SYNPD1 vs SYNPD2
        group.name2=models(contains(models,'PD'));
        
        model1=Allmodels.(type{i}).(group.name1{1}); % SYN1
        model2=Allmodels.(type{i}).(group.name1{2}); % SYN2
        model3=Allmodels.(type{i}).(group.name2{1}); % SYNPD1
        model4=Allmodels.(type{i}).(group.name2{2}); % SYNPD2
        
        % add mitochondrial rxns/genes
        model1.mrxns=findRxnFromCompartment(model1,'[m]');
        model2.mrxns=findRxnFromCompartment(model2,'[m]');
        model3.mrxns=findRxnFromCompartment(model3,'[m]');
        model4.mrxns=findRxnFromCompartment(model4,'[m]');
        
        subnetwork=extractSubNetwork(model1,model1.mrxns);
        subnetwork=updateGenes(subnetwork);
        model1.mgenes=subnetwork.genes;
        subnetwork=extractSubNetwork(model2,model2.mrxns);
        subnetwork=updateGenes(subnetwork);
        model2.mgenes=subnetwork.genes;
        subnetwork=extractSubNetwork(model3,model3.mrxns);
        subnetwork=updateGenes(subnetwork);
        model3.mgenes=subnetwork.genes;
        subnetwork=extractSubNetwork(model4,model4.mrxns);
        subnetwork=updateGenes(subnetwork);
        model4.mgenes=subnetwork.genes;

        
        %%%%%%%%%%%%%% total all rxns
        compareExRxnsTable.(type{i}).TotalRxns=cell(5,4);
        compareExRxnsTable.(type{i}).TotalRxns{1,1}='crossmatch on total rxns';
        compareExRxnsTable.(type{i}).TotalRxns{1,2}='model1';
        compareExRxnsTable.(type{i}).TotalRxns{1,3}='model2';
        compareExRxnsTable.(type{i}).TotalRxns{1,4}='intersect';
        % SYN1 vs SYN2
        compareExRxnsTable.(type{i}).TotalRxns{2,1}=[group.name1{1} ' vs ' group.name1{2}];
        compareExRxnsTable.(type{i}).TotalRxns{2,2}=length(model1.rxns);
        compareExRxnsTable.(type{i}).TotalRxns{2,3}=length(model2.rxns);
        compareExRxnsTable.(type{i}).TotalRxns{2,4}=length(intersect(model1.rxns,model2.rxns));
        % SYNPD1 vs SYNPD2
        compareExRxnsTable.(type{i}).TotalRxns{3,1}=[group.name2{1} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).TotalRxns{3,2}=length(model3.rxns);
        compareExRxnsTable.(type{i}).TotalRxns{3,3}=length(model4.rxns);
        compareExRxnsTable.(type{i}).TotalRxns{3,4}=length(intersect(model3.rxns,model4.rxns));
        % SYN1 vs SYNPD1
        compareExRxnsTable.(type{i}).TotalRxns{4,1}=[group.name1{1} ' vs ' group.name2{1}];
        compareExRxnsTable.(type{i}).TotalRxns{4,2}=length(model1.rxns);
        compareExRxnsTable.(type{i}).TotalRxns{4,3}=length(model3.rxns);
        compareExRxnsTable.(type{i}).TotalRxns{4,4}=length(intersect(model1.rxns,model3.rxns));
        % SYN2 vs SYNPD2
        compareExRxnsTable.(type{i}).TotalRxns{5,1}=[group.name1{2} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).TotalRxns{5,2}=length(model2.rxns);
        compareExRxnsTable.(type{i}).TotalRxns{5,3}=length(model4.rxns);
        compareExRxnsTable.(type{i}).TotalRxns{5,4}=length(intersect(model2.rxns,model4.rxns));
        
        % add bool        
        model1 = findSExRxnInd(model1);
        model2 = findSExRxnInd(model2);
        model3 = findSExRxnInd(model3);
        model4 = findSExRxnInd(model4);
           
        [model1.synExc,model1.synUpt] = findExcRxns(model1);
        [model2.synExc,model2.synUpt] = findExcRxns(model2);
        [model3.synExc,model3.synUpt] = findExcRxns(model3);
        [model4.synExc,model4.synUpt] = findExcRxns(model4);
        
        %%%%%%%%%%%%%% total Exchange rxns
        compareExRxnsTable.(type{i}).ExRxns=cell(5,4);
        compareExRxnsTable.(type{i}).ExRxns{1,1}='crossmatch for Exchange/SInk/DM rxns';
        compareExRxnsTable.(type{i}).ExRxns{1,2}='model1';
        compareExRxnsTable.(type{i}).ExRxns{1,3}='model2';
        compareExRxnsTable.(type{i}).ExRxns{1,4}='intersect';
        % SYN1 vs SYN2
        compareExRxnsTable.(type{i}).ExRxns{2,1}=[group.name1{1} ' vs ' group.name1{2}];
        compareExRxnsTable.(type{i}).ExRxns{2,2}=sum(model1.ExchRxnBool | model1.SinkRxnBool | model1.DMRxnBool);
        compareExRxnsTable.(type{i}).ExRxns{2,3}=sum(model2.ExchRxnBool | model2.SinkRxnBool | model2.DMRxnBool);
        compareExRxnsTable.(type{i}).ExRxns{2,4}=length(intersect(model1.rxns(model1.ExchRxnBool | model1.SinkRxnBool | model1.DMRxnBool),model2.rxns(model2.ExchRxnBool | model2.SinkRxnBool | model2.DMRxnBool)));
        % SYNPD1 vs SYNPD2
        compareExRxnsTable.(type{i}).ExRxns{3,1}=[group.name2{1} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).ExRxns{3,2}=sum(model3.ExchRxnBool | model3.SinkRxnBool | model3.DMRxnBool);
        compareExRxnsTable.(type{i}).ExRxns{3,3}=sum(model4.ExchRxnBool | model4.SinkRxnBool | model4.DMRxnBool);
        compareExRxnsTable.(type{i}).ExRxns{3,4}=length(intersect(model3.rxns(model3.ExchRxnBool | model3.SinkRxnBool | model3.DMRxnBool),model4.rxns(model4.ExchRxnBool | model4.SinkRxnBool | model4.DMRxnBool)));
        % SYN1 vs SYNPD1
        compareExRxnsTable.(type{i}).ExRxns{4,1}=[group.name1{1} ' vs ' group.name2{1}];
        compareExRxnsTable.(type{i}).ExRxns{4,2}=sum(model1.ExchRxnBool | model1.SinkRxnBool | model1.DMRxnBool);
        compareExRxnsTable.(type{i}).ExRxns{4,3}=sum(model3.ExchRxnBool | model3.SinkRxnBool | model3.DMRxnBool);
        compareExRxnsTable.(type{i}).ExRxns{4,4}=length(intersect(model1.rxns(model1.ExchRxnBool | model1.SinkRxnBool | model1.DMRxnBool),model3.rxns(model3.ExchRxnBool | model3.SinkRxnBool | model3.DMRxnBool)));
        % SYN2 vs SYNPD2
        compareExRxnsTable.(type{i}).ExRxns{5,1}=[group.name1{2} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).ExRxns{5,2}=sum(model2.ExchRxnBool | model2.SinkRxnBool | model2.DMRxnBool);
        compareExRxnsTable.(type{i}).ExRxns{5,3}=sum(model4.ExchRxnBool | model4.SinkRxnBool | model4.DMRxnBool);
        compareExRxnsTable.(type{i}).ExRxns{5,4}=length(intersect(model2.rxns(model2.ExchRxnBool | model2.SinkRxnBool | model2.DMRxnBool),model4.rxns(model4.ExchRxnBool | model4.SinkRxnBool | model4.DMRxnBool)));
        
        %%%%%%%%%%%%%% total mitochondrial rxns
        compareExRxnsTable.(type{i}).mitoRxns=cell(5,4);
        compareExRxnsTable.(type{i}).mitoRxns{1,1}='crossmatch on mitochondrial rxns';
        compareExRxnsTable.(type{i}).mitoRxns{1,2}='model1';
        compareExRxnsTable.(type{i}).mitoRxns{1,3}='model2';
        compareExRxnsTable.(type{i}).mitoRxns{1,4}='intersect';
        % SYN1 vs SYN2
        compareExRxnsTable.(type{i}).mitoRxns{2,1}=[group.name1{1} ' vs ' group.name1{2}];
        compareExRxnsTable.(type{i}).mitoRxns{2,2}=length(model1.mrxns);
        compareExRxnsTable.(type{i}).mitoRxns{2,3}=length(model2.mrxns);
        compareExRxnsTable.(type{i}).mitoRxns{2,4}=length(intersect(model1.mrxns(:,1),model2.mrxns(:,1)));
        % SYNPD1 vs SYNPD2
        compareExRxnsTable.(type{i}).mitoRxns{3,1}=[group.name2{1} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).mitoRxns{3,2}=length(model3.mrxns);
        compareExRxnsTable.(type{i}).mitoRxns{3,3}=length(model4.mrxns);
        compareExRxnsTable.(type{i}).mitoRxns{3,4}=length(intersect(model3.mrxns(:,1),model4.mrxns(:,1)));
        % SYN1 vs SYNPD1
        compareExRxnsTable.(type{i}).mitoRxns{4,1}=[group.name1{1} ' vs ' group.name2{1}];
        compareExRxnsTable.(type{i}).mitoRxns{4,2}=length(model1.mrxns);
        compareExRxnsTable.(type{i}).mitoRxns{4,3}=length(model3.mrxns);
        compareExRxnsTable.(type{i}).mitoRxns{4,4}=length(intersect(model1.mrxns(:,1),model3.mrxns(:,1)));
        % SYN2 vs SYNPD2
        compareExRxnsTable.(type{i}).mitoRxns{5,1}=[group.name1{2} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).mitoRxns{5,2}=length(model2.mrxns);
        compareExRxnsTable.(type{i}).mitoRxns{5,3}=length(model4.mrxns);
        compareExRxnsTable.(type{i}).mitoRxns{5,4}=length(intersect(model2.mrxns(:,1),model4.mrxns(:,1)));
        
        %%%%%%%%%%%%%% total genes
        compareExRxnsTable.(type{i}).Totalgenes=cell(5,4);
        compareExRxnsTable.(type{i}).Totalgenes{1,1}='crossmatch on total genes';
        compareExRxnsTable.(type{i}).Totalgenes{1,2}='model1';
        compareExRxnsTable.(type{i}).Totalgenes{1,3}='model2';
        compareExRxnsTable.(type{i}).Totalgenes{1,4}='intersect';
        % SYN1 vs SYN2
        compareExRxnsTable.(type{i}).Totalgenes{2,1}=[group.name1{1} ' vs ' group.name1{2}];
        compareExRxnsTable.(type{i}).Totalgenes{2,2}=length(model1.genes);
        compareExRxnsTable.(type{i}).Totalgenes{2,3}=length(model2.genes);
        compareExRxnsTable.(type{i}).Totalgenes{2,4}=length(intersect(model1.genes,model2.genes));
        % SYNPD1 vs SYNPD2
        compareExRxnsTable.(type{i}).Totalgenes{3,1}=[group.name2{1} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).Totalgenes{3,2}=length(model3.genes);
        compareExRxnsTable.(type{i}).Totalgenes{3,3}=length(model4.genes);
        compareExRxnsTable.(type{i}).Totalgenes{3,4}=length(intersect(model3.genes,model4.genes));
        % SYN1 vs SYNPD1
        compareExRxnsTable.(type{i}).Totalgenes{4,1}=[group.name1{1} ' vs ' group.name2{1}];
        compareExRxnsTable.(type{i}).Totalgenes{4,2}=length(model1.genes);
        compareExRxnsTable.(type{i}).Totalgenes{4,3}=length(model3.genes);
        compareExRxnsTable.(type{i}).Totalgenes{4,4}=length(intersect(model1.genes,model3.genes));
        % SYN2 vs SYNPD2
        compareExRxnsTable.(type{i}).Totalgenes{5,1}=[group.name1{2} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).Totalgenes{5,2}=length(model2.genes);
        compareExRxnsTable.(type{i}).Totalgenes{5,3}=length(model4.genes);
        compareExRxnsTable.(type{i}).Totalgenes{5,4}=length(intersect(model2.genes,model4.genes));
        
        %%%%%%%%%%%%%% total mitochondrial genes
        compareExRxnsTable.(type{i}).mitogenes=cell(5,4);
        compareExRxnsTable.(type{i}).mitogenes{1,1}='crossmatch on total genes';
        compareExRxnsTable.(type{i}).mitogenes{1,2}='model1';
        compareExRxnsTable.(type{i}).mitogenes{1,3}='model2';
        compareExRxnsTable.(type{i}).mitogenes{1,4}='intersect';
        % SYN1 vs SYN2
        compareExRxnsTable.(type{i}).mitogenes{2,1}=[group.name1{1} ' vs ' group.name1{2}];
        compareExRxnsTable.(type{i}).mitogenes{2,2}=length(model1.mgenes);
        compareExRxnsTable.(type{i}).mitogenes{2,3}=length(model2.mgenes);
        compareExRxnsTable.(type{i}).mitogenes{2,4}=length(intersect(model1.mgenes,model2.mgenes));
        % SYNPD1 vs SYNPD2
        compareExRxnsTable.(type{i}).mitogenes{3,1}=[group.name2{1} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).mitogenes{3,2}=length(model3.mgenes);
        compareExRxnsTable.(type{i}).mitogenes{3,3}=length(model4.mgenes);
        compareExRxnsTable.(type{i}).mitogenes{3,4}=length(intersect(model3.mgenes,model4.mgenes));
        % SYN1 vs SYNPD1
        compareExRxnsTable.(type{i}).mitogenes{4,1}=[group.name1{1} ' vs ' group.name2{1}];
        compareExRxnsTable.(type{i}).mitogenes{4,2}=length(model1.mgenes);
        compareExRxnsTable.(type{i}).mitogenes{4,3}=length(model3.mgenes);
        compareExRxnsTable.(type{i}).mitogenes{4,4}=length(intersect(model1.mgenes,model3.mgenes));
        % SYN2 vs SYNPD2
        compareExRxnsTable.(type{i}).mitogenes{5,1}=[group.name1{2} ' vs ' group.name2{2}];
        compareExRxnsTable.(type{i}).mitogenes{5,2}=length(model2.mgenes);
        compareExRxnsTable.(type{i}).mitogenes{5,3}=length(model4.mgenes);
        compareExRxnsTable.(type{i}).mitogenes{5,4}=length(intersect(model2.mgenes,model4.mgenes));
end
end

