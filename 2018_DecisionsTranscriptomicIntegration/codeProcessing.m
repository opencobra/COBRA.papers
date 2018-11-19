%This code is related to the following publication
%"Assessing key decisions for transcriptomic data integration in
%biochemical networks"
%Authors : Anne Richelle, Chintan Joshi and Nathan E Lewis
%doi: https://doi.org/10.1101/301945
%Please cite it if you use it for your research
%If you have any question, do not hesitate to contact me:
%arichelleres@gmail.com

%%____________________________________________________________________________________________________________%%
%------------ANALYSIS OF TRANSCRIPTOMICS DATA AND RECON2.2 MODEL CONTENT---
%%____________________________________________________________________________________________________________%%

% download the model Recon 2.2 from https://www.ebi.ac.uk/biomodels/MODEL1603150001
system('curl -O https://www.ebi.ac.uk/biomodels/model/download/MODEL1603150001.2\?filename\=MODEL1603150001_url.xml')
system('mv MODEL1603150001.2\?filename\=MODEL1603150001_url.xml Recon2.2.xml')

%% Identify the metabolic genes of HPA present in recon2.2 model
metabolicGene=[];
readCbModel('Recon2.2.xml')
for i=1:length(model.genes)
    metabolicGene(i)=str2num(model.genes{i});
end
metabolicGene=unique(metabolicGene);

load('dataHPA')
data_extraction=[];
missing=[];
for i=1:length(metabolicGene)
    ID=find(data(:,1)==metabolicGene(i));
    if isempty(ID)
        missing(end+1)=metabolicGene(i);
    else
        data_extraction(end+1,:)=data(ID,:);
    end
end

data = data_extraction;
expressionData={};
expressionData.gene=data(:,1);
expressionData.mean_value=data(:,2);
expressionData.valuebyTissue=data(:,3:end);
expressionData.min=min(expressionData.valuebyTissue,[],2);
expressionData.Tissue={'adipose tissue','adrenal gland','appendix','bone marrow','brain','colon','duodenum','endometrium','esophagus','fallopian tube','gallbladder','heart muscle','kidney','liver','lung','lymph node','ovary','pancreas','placenta','prostate','rectum','salivary gland','skeletal muscle','skin','small intestine','smooth muscle','spleen','stomach','testis','thyroid gland','tonsil','urinary bladder'};


% Removal of the missing genes from the model
ID_geneMissing=[];
for i=1:length(missing)
    ID_geneMissing(i)=findGeneIDs(model,num2str(missing(i)));
end
model= removeFieldEntriesForType(model,ID_geneMissing,'genes', 1675);
model = creategrRulesField(model);

%% - Distribution of gene expression value for HPA dataset mapped to recon2.2
% compute the distribution of gene expression observed in global data (only metabolic gene) to
% define a global threshold of gene activity

expressionGlobal=[];
for i=1:size(expressionData.valuebyTissue,2)
    expressionGlobal=[expressionGlobal;expressionData.valuebyTissue(:,i)];
end

expressionGlobal(expressionGlobal<=-max(expressionGlobal)) = -max(expressionGlobal);
exprG=log10(expressionGlobal);

figure(1);hist(exprG,50)
    l95= (prctile(exprG,95));
    ths_95=10^l95;
    l90= (prctile(exprG,90));
    ths_90=10^l90;
    l75 = (prctile(exprG,75));
    ths_75=10^l75;
    l50 = (prctile(exprG,50));
    ths_50=10^l50;
    l25 = (prctile(exprG,25));
    ths_25=10^l25;
    l10 = (prctile(exprG,10));
    ths_10=10^l10;
    h95=line([l95,l95],get(gca,'YLim'),'Color','m');
    h90=line([l90,l90],get(gca,'YLim'),'Color','y');
    h75=line([l75,l75],get(gca,'YLim'),'Color','b');
    h50=line([l50,l50],get(gca,'YLim'),'Color','c');
    h25=line([l25,l25],get(gca,'YLim'),'Color','g');
    h10=line([l10,l10],get(gca,'YLim'),'Color','r');
legend([h95,h90,h75,h50,h25,h10],{['95th percentile = ',num2str(ths_95),' FPKM'],['90th percentile = ',num2str(ths_90),' FPKM'],['75th percentile = ',num2str(ths_75),' FPKM'],['50th percentile = ',num2str(ths_50),' FPKM'],['25th percentile = ',num2str(ths_25),' FPKM'],['10th percentile = ',num2str(ths_10),' FPKM']},'FontSize',14)
xlabel('log10(expressionValue)','FontSize',14)
ylabel('Genes','FontSize',14)

%% Computation of the threshold values for T1 and T2 states definition in local and global approaches
%_______GLOBAL_____%
expression.ths_95=10^l95;
expression.ths_90=10^l90;
expression.ths_75=10^l75;
expression.ths_50=10^l50;
expression.ths_25=10^l25;
expression.ths_10=10^l10;

%_______LOCAL_____%
%------ Local T1 - 25th------%
expression.ths_local_T1_25=[];
for i=1:length(expressionData.gene)
    	expressionValue=expressionData.valuebyTissue(i,:);
      	expression.ths_local_T1_25(i,:)=max(mean(expressionValue),expression.ths_25);
end
%------ Local T2 - 25th & 75th------%
expression.ths_local_T2_25_75=[];
for i=1:length(expressionData.gene)
    	expressionValue=expressionData.valuebyTissue(i,:);
        if  mean(expressionValue)>= expression.ths_75
            expression.ths_local_T2_25_75(i)=expression.ths_75;
        else
            expression.ths_local_T2_25_75(i)=max(mean(expressionValue),expression.ths_25);
        end
end
%------ Local T2 - 25th & 90th------%
expression.ths_local_T2_25_90=[];
for i=1:length(expressionData.gene)
    	expressionValue=expressionData.valuebyTissue(i,:);
        if  mean(expressionValue)>= expression.ths_90
            expression.ths_local_T2_25_90(i)=expression.ths_90;
        else
            expression.ths_local_T2_25_90(i)=max(mean(expressionValue),expression.ths_25);
        end
end

%% Combinations of preprocessing decisions
expression.scoreGlobal=expressionData.valuebyTissue;
expression.Tissue=expressionData.Tissue;

%_________________________________________________________________________________________________________%
%_____________________________________________C1. CASE1___________________________________________________%
%_________________________________________________________________________________________________________%
% Perform first the thresholding on all the metabolic gene and after the
% gene mapping
%_____________________________________________C1.1 - GLOBAL________________________________________________%
%------------------------------------------C1.1.1 -  T1:50th-----------------------------------------------%
%--------------------------------------------C1.1.1.1 - GM1------------------------------------------------%
%Note: GM1 = An AND will be replaced by MIN and an OR will be replaced by MAX.
display('Case 2 - GM1 - Global T1 50th')
scoreUp50=expression.scoreGlobal;
scoreUp50(scoreUp50<=expression.ths_50)=0;
expression.Rxns_case2_GM1_global50=[];
expression.geneUsed_case2_GM1_global50={};
for i=1:length(expressionData.Tissue)
    expressionData.value=scoreUp50(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'false');
    expression.Rxns_case2_GM1_global50=[expression.Rxns_case2_GM1_global50 expressionRxns];
    expression.geneUsed_case2_GM1_global50{i}= gene_used;
end
%--------------------------------------------C1.1.1.2 - GM2------------------------------------------------%
%Note: GM2 = An AND will be replaced by MIN and an OR will be replaced by SUM.
display('Case 2 - GM2 - Global T1 50th')
expression.Rxns_case2_GM2_global50=[];
expression.geneUsed_case2_GM2_global50={};
for i=1:length(expressionData.Tissue)
    expressionData.value=scoreUp50(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'true');
    expression.Rxns_case2_GM2_global50=[expression.Rxns_case2_GM2_global50 expressionRxns];
    expression.geneUsed_case2_GM2_global50{i}= gene_used;
end
%------------------------------------------C1.1.2 - T1:75th-----------------------------------------------%
%-------------------------------------------C1.1.2.1 - GM1------------------------------------------------%
%Note: GM1 = An AND will be replaced by MIN and an OR will be replaced by MAX.
display('Case 2 - GM1 - Global T1 75th')
scoreUp75=expression.scoreGlobal;
scoreUp75(scoreUp75<=expression.ths_75)=0;
expression.Rxns_case2_GM1_global75=[];
expression.geneUsed_case2_GM1_global75={};
for i=1:length(expressionData.Tissue)
    expressionData.value=scoreUp75(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'false');
    expression.Rxns_case2_GM1_global75=[expression.Rxns_case2_GM1_global75 expressionRxns];
    expression.geneUsed_case2_GM1_global75{i}= gene_used;
end
%-------------------------------------------C1.1.2.2 - GM2------------------------------------------------%
%Note: GM2 = An AND will be replaced by MIN and an OR will be replaced by SUM.
display('Case 2 - GM2 - Global T1 75th')
expression.Rxns_case2_GM2_global75=[];
expression.geneUsed_case2_GM2_global75={};
for i=1:length(expressionData.Tissue)
    expressionData.value=scoreUp75(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'true');
    expression.Rxns_case2_GM2_global75=[expression.Rxns_case2_GM2_global75 expressionRxns];
    expression.geneUsed_case2_GM2_global75{i}= gene_used;
end
%_____________________________________________C1.2 - LOCAL________________________________________________%
%-------------------------------------------C1.2.1 - T1:25th----------------------------------------------%
expression.scoreLocal_T1_25=[];
for i=1:length(expressionData.Tissue)
    expression.scoreLocal_T1_25(:,i)=(expression.scoreGlobal(:,i)./expression.ths_local_T1_25)';
end
%-------------------------------------------C1.2.1.1 - GM1----------------------------------------------%
%Note: GM1 = An AND will be replaced by MIN and an OR will be replaced by MAX.
display('Case 2 - GM1 - Local T1 25th')
expression.Rxns_case2_GM1_local25=[];
expression.geneUsed_case2_GM1_local25={};
for i=1:length(expressionData.Tissue)
    expressionData.value=expression.scoreLocal_T1_25(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'false');
    expression.Rxns_case2_GM1_local25=[expression.Rxns_case2_GM1_local25 expressionRxns];
    expression.geneUsed_case2_GM1_local25{i}= gene_used;
end
%-------------------------------------------C1.2.1.2 - GM2----------------------------------------------%
%Note: GM2 = An AND will be replaced by MIN and an OR will be replaced by SUM.
display('Case 2 - GM2 - Local T1 25th')
expression.Rxns_case2_GM2_local25=[];
expression.geneUsed_case2_GM2_local25={};
for i=1:length(expressionData.Tissue)
    expressionData.value=expression.scoreLocal_T1_25(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'true');
    expression.Rxns_case2_GM2_local25=[expression.Rxns_case2_GM2_local25 expressionRxns];
    expression.geneUsed_case2_GM2_local25{i}= gene_used;
end
%-------------------------------------------C1.2.2 - T2:25th & 75th-------------------------------------------%
expression.scoreLocal_T2_25_75=[];
for i=1:length(expressionData.Tissue)
    expression.scoreLocal_T2_25_75(:,i)=(expression.scoreGlobal(:,i)./expression.ths_local_T2_25_75');
end
%-------------------------------------------C1.2.2.1 - GM1--------------------------------------------------%
%Note: GM1 = An AND will be replaced by MIN and an OR will be replaced by MAX.
display('Case 2 - GM1 - Local 25th & 75th')
expression.Rxns_case2_GM1_local25_75=[];
expression.geneUsed_case2_GM1_local25_75={};
for i=1:length(expressionData.Tissue)
    expressionData.value=expression.scoreLocal_T2_25_75(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'false');
    expression.Rxns_case2_GM1_local25_75=[expression.Rxns_case2_GM1_local25_75 expressionRxns];
    expression.geneUsed_case2_GM1_local25_75{i}= gene_used;
end
%-------------------------------------------C1.2.2.2 - GM2-------------------------------------------------%
%Note: GM2 = An AND will be replaced by MIN and an OR will be replaced by SUM.
display('Case 2 - GM2 - Local 25th & 75th')
expression.Rxns_case2_GM2_local25_75=[];
expression.geneUsed_case2_GM2_local25_75={};
for i=1:length(expressionData.Tissue)
    expressionData.value=expression.scoreLocal_T2_25_75(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'true');
    expression.Rxns_case2_GM2_local25_75=[expression.Rxns_case2_GM2_local25_75 expressionRxns];
    expression.geneUsed_case2_GM2_local25_75{i}= gene_used;
end
%-------------------------------------------C1.2.3 - T2:25th & 90th-------------------------------------------%
expression.scoreLocal_T2_25_90=[];
for i=1:length(expressionData.Tissue)
    expression.scoreLocal_T2_25_90(:,i)=(expression.scoreGlobal(:,i)./expression.ths_local_T2_25_90');
end
%-------------------------------------------C1.2.3.1 - GM1-------------------------------------------%
%Note: GM1 = An AND will be replaced by MIN and an OR will be replaced by MAX.
display('Case 2 - GM1 - Local 25th & 90th')
expression.Rxns_case2_GM1_local25_90=[];
expression.geneUsed_case2_GM1_local25_90={};
for i=1:length(expressionData.Tissue)
    expressionData.value=expression.scoreLocal_T2_25_90(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'false');
    expression.Rxns_case2_GM1_local25_90=[expression.Rxns_case2_GM1_local25_90 expressionRxns];
    expression.geneUsed_case2_GM1_local25_90{i}= gene_used;
end
%-------------------------------------------C1.2.3.2 - GM2-------------------------------------------%
%Note: GM2 = An AND will be replaced by MIN and an OR will be replaced by SUM.
display('Case 2 - GM2 - Local 25th & 90th')
expression.Rxns_case2_GM2_local25_90=[];
expression.geneUsed_case2_GM2_local25_90={};
for i=1:length(expressionData.Tissue)
    expressionData.value=expression.scoreLocal_T2_25_90(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'true');
    expression.Rxns_case2_GM2_local25_90=[expression.Rxns_case2_GM2_local25_90 expressionRxns];
    expression.geneUsed_case2_GM2_local25_90{i}= gene_used;
end
%_________________________________________________________________________________________________________%
%_____________________________________________C2. CASE1___________________________________________________%
%_________________________________________________________________________________________________________%
% Perform first the gene mapping and after the thresholding on mapped genes
%---------------------------------------------C2.1 - GM1-------------------------------------------------%
%Note: GM1 = An AND will be replaced by MIN and an OR will be replaced by MAX.
display('Case 1 - GM1')
expression.Rxns_case1_GM1=[];
expression.geneUsed_case1_GM1={};
for i=1:length(expressionData.Tissue)
    expressionData.value=expression.scoreGlobal(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'false');
    expression.Rxns_case1_GM1=[expression.Rxns_case1_GM1 expressionRxns];
    expression.geneUsed_case1_GM1{i}= gene_used;
end
expression.ths_case1_GM1_local25=[];
expression.ths_case1_GM1_local25_75=[];
expression.ths_case1_GM1_local25_90=[];
for k=1:length(expressionData.Tissue)
    geneUsed_case1_GM1=expression.geneUsed_case1_GM1{k};
    for i=1:7785
        if ~isempty(geneUsed_case1_GM1{i})
            gene=geneUsed_case1_GM1{i};
            ID_case1_GM1=find(expressionData.gene==str2num(gene{1}));
            expressionValue=expressionData.valuebyTissue(ID_case1_GM1,:);
%---------------------------------------------C2.1.1 - Local T1:25th------------------------------------%
            expression.ths_case1_GM1_local25(i,k)=max(mean(expressionValue),expression.ths_25);
%-------------------------------------------C2.1.2 - Local T1:25th & 75th------------------------------%
            if  mean(expressionValue)>= expression.ths_75
                expression.ths_case1_GM1_local25_75(i,k)=expression.ths_75;
            else
                expression.ths_case1_GM1_local25_75(i,k)=max(mean(expressionValue),expression.ths_25);
            end
%-------------------------------------------C2.1.3 - Local T1:25th & 90th------------------------------%
            if  mean(expressionValue)>= expression.ths_90
                expression.ths_case1_GM1_local25_90(i,k)=expression.ths_90;
            else
                expression.ths_case1_GM1_local25_90(i,k)=max(mean(expressionValue),expression.ths_25);
            end
        else
            expression.ths_case1_GM1_local25(i,k)=1;
            expression.ths_case1_GM1_local25_75(i,k)=1;
            expression.ths_case1_GM1_local25_90(i,k)=1;
        end
    end
end
%---------------------------------------------C2.2 - GM2-------------------------------------------------%
%Note: GM2 = An AND will be replaced by MIN and an OR will be replaced by SUM.
display('Case 1 - GM2')
expression.Rxns_case1_GM2=[];
expression.geneUsed_case1_GM2={};
for i=1:length(expressionData.Tissue)
    expressionData.value=expression.scoreGlobal(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model,expressionData,'true');
    expression.Rxns_case1_GM2= [expression.Rxns_case1_GM2 expressionRxns];
    expression.geneUsed_case1_GM2= gene_used;
end
expression.ths_case1_GM2_local25=[];
expression.ths_case1_GM2_local25_75=[];
expression.ths_case1_GM2_local25_90=[];
geneUsed_case1_GM2=expression.geneUsed_case1_GM2;
for i=1:7785
	if ~isempty(geneUsed_case1_GM2{i})
        gene=geneUsed_case1_GM2{i};
        local25=[];
        local25_75=[];
        local25_90=[];
        for j=1:length(geneUsed_case1_GM2{i})
            ID=find(expressionData.gene==str2num(gene{j}));
            expressionValue=expressionData.valuebyTissue(ID,:);
%---------------------------------------------C2.2.1 - Local T1:25th------------------------------------%
            local25(j)=max(mean(expressionValue),expression.ths_25);
%-----------------------------------------C2.2.2 - Local T2:25th & 75th------------------------------------%
            if  mean(expressionValue)>= expression.ths_75
                local25_75(j)=expression.ths_75;
            else
                local25_75(j)=max(mean(expressionValue),expression.ths_25);
            end
%----------------------------------------C2.2.3 - Local T2:25th & 90th------------------------------------%
            if  mean(expressionValue)>= expression.ths_90
                local25_90(j)=expression.ths_90;
            else
                local25_90(j)=max(mean(expressionValue),expression.ths_25);
            end
        end
        expression.ths_case1_GM2_local25(i)=sum(local25);
        expression.ths_case1_GM2_local25_75(i)=sum(local25_75);
        expression.ths_case1_GM2_local25_90(i)=sum(local25_90);
    else
     	expression.ths_case1_GM2_local25(i)=1;
     	expression.ths_case1_GM2_local25_75(i)=1;
    	expression.ths_case1_GM2_local25_90(i)=1;
	end
end

%% How many reactions are considered as active depending on preprocessing combinations used
countActiveRxns_case1_GM1_global50=[];
countActiveRxns_case1_GM2_global50=[];
countActiveRxns_case1_GM1_global75=[];
countActiveRxns_case1_GM2_global75=[];
countActiveRxns_case1_GM1_local25=[];
countActiveRxns_case1_GM2_local25=[];
countActiveRxns_case1_GM1_local25_75=[];
countActiveRxns_case1_GM2_local25_75=[];
countActiveRxns_case1_GM1_local25_90=[];
countActiveRxns_case1_GM2_local25_90=[];
countActiveRxns_case2_GM1_global50=[];
countActiveRxns_case2_GM2_global50=[];
countActiveRxns_case2_GM1_global75=[];
countActiveRxns_case2_GM2_global75=[];
countActiveRxns_case2_GM1_local25=[];
countActiveRxns_case2_GM2_local25=[];
countActiveRxns_case2_GM1_local25_75=[];
countActiveRxns_case2_GM2_local25_75=[];
countActiveRxns_case2_GM1_local25_90=[];
countActiveRxns_case2_GM2_local25_90=[];

for i=1:length(expressionData.Tissue)
    countActiveRxns_case1_GM1_global50(i)=length(find(expression.Rxns_case1_GM1(:,i)>=expression.ths_50));
    countActiveRxns_case1_GM2_global50(i)=length(find(expression.Rxns_case1_GM2(:,i)>=expression.ths_50));
    countActiveRxns_case2_GM1_global50(i)=length(find(expression.Rxns_case2_GM1_global50(:,i)>=expression.ths_50));
    countActiveRxns_case2_GM2_global50(i)=length(find(expression.Rxns_case2_GM2_global50(:,i)>=expression.ths_50));
    countActiveRxns_case1_GM1_global75(i)=length(find(expression.Rxns_case1_GM1(:,i)>=expression.ths_75));
    countActiveRxns_case1_GM2_global75(i)=length(find(expression.Rxns_case1_GM2(:,i)>=expression.ths_75));
    countActiveRxns_case2_GM1_global75(i)=length(find(expression.Rxns_case2_GM1_global75(:,i)>=expression.ths_75));
    countActiveRxns_case2_GM2_global75(i)=length(find(expression.Rxns_case2_GM2_global75(:,i)>=expression.ths_75));
    countActiveRxns_case1_GM1_local25(i)=length(find((expression.Rxns_case1_GM1(:,i)./expression.ths_case1_GM1_local25(:,i))>=1));
    countActiveRxns_case1_GM2_local25(i)=length(find((expression.Rxns_case1_GM2(:,i)./expression.ths_case1_GM2_local25')>=1));
    countActiveRxns_case1_GM1_local25_75(i)=length(find((expression.Rxns_case1_GM1(:,i)./expression.ths_case1_GM1_local25_75(:,i))>=1));
    countActiveRxns_case1_GM2_local25_75(i)=length(find((expression.Rxns_case1_GM2(:,i)./expression.ths_case1_GM2_local25_75')>=1));
    countActiveRxns_case1_GM1_local25_90(i)=length(find((expression.Rxns_case1_GM1(:,i)./expression.ths_case1_GM1_local25_90(:,i))>=1));
    countActiveRxns_case1_GM2_local25_90(i)=length(find((expression.Rxns_case1_GM2(:,i)./expression.ths_case1_GM2_local25_90')>=1));
    countActiveRxns_case2_GM1_local25(i)=length(find(expression.Rxns_case2_GM1_local25(:,i)>=1));
    countActiveRxns_case2_GM2_local25(i)=length(find(expression.Rxns_case2_GM2_local25(:,i)>=1));
    countActiveRxns_case2_GM1_local25_75(i)=length(find(expression.Rxns_case2_GM1_local25_75(:,i)>=1));
    countActiveRxns_case2_GM2_local25_75(i)=length(find(expression.Rxns_case2_GM2_local25_75(:,i)>=1));
    countActiveRxns_case2_GM1_local25_90(i)=length(find(expression.Rxns_case2_GM1_local25_90(:,i)>=1));
    countActiveRxns_case2_GM2_local25_90(i)=length(find(expression.Rxns_case2_GM2_local25_90(:,i)>=1));
end

ActiveRxnscount=[countActiveRxns_case1_GM1_global50' ...
countActiveRxns_case2_GM1_global50' ...
countActiveRxns_case1_GM2_global50' ...
countActiveRxns_case2_GM2_global50' ...
countActiveRxns_case1_GM1_global75' ...
countActiveRxns_case2_GM1_global75' ...
countActiveRxns_case1_GM2_global75' ...
countActiveRxns_case2_GM2_global75' ...
countActiveRxns_case1_GM1_local25' ...
countActiveRxns_case2_GM1_local25' ...
countActiveRxns_case1_GM2_local25' ...
countActiveRxns_case2_GM2_local25' ...
countActiveRxns_case1_GM1_local25_75' ...
countActiveRxns_case2_GM1_local25_75' ...
countActiveRxns_case1_GM2_local25_75' ...
countActiveRxns_case2_GM2_local25_75' ...
countActiveRxns_case1_GM1_local25_90' ...
countActiveRxns_case2_GM1_local25_90' ...
countActiveRxns_case1_GM2_local25_90' ...
countActiveRxns_case2_GM2_local25_90'];

count_percent=ActiveRxnscount./4724*100;

figure(2);
boxplot(count_percent)
ax=gca;
ax.XTickLabel={'Case 1 GM1 Global T1 50th' ...
'Case 2 GM1 Global T1 50th' ...
'Case 1 GM2 Global T1 50th' ...
'Case 2 GM2 Global T1 50th' ...
'Case 1 GM1 Global T1 75th' ...
'Case 2 GM1 Global T1 75th' ...
'Case 1 GM2 Global T1 75th' ...
'Case 2 GM2 Global T1 75th' ...
'Case 1 GM1 Local T1 25th' ...
'Case 2 GM1 Local T1 25th' ...
'Case 1 GM2 Local T1 25th' ...
'Case 2 GM2 Local T1 25th' ...
'Case 1 GM1 Local T2 25th & 75th' ...
'Case 2 GM1 Local T2 25th & 75th' ...
'Case 1 GM2 Local T2 25th & 75th' ...
'Case 2 GM2 Local T2 25th & 75th' ...
'Case 1 GM1 Local T2 25th & 90th' ...
'Case 2 GM1 Local T2 25th & 90th' ...
'Case 1 GM2 Local T2 25th & 90th' ...
'Case 2 GM2 Local T2 25th & 90th'};
ax.XTickLabelRotation=45;
ylabel('Percentage of reactions considered as active','FontSize',14)

%% What is the similarity between reactions considered as active depending on preprocessing combinations
ActiveRxns_case1_GM1_global50={};
ActiveRxns_case1_GM2_global50={};
ActiveRxns_case1_GM1_global75={};
ActiveRxns_case1_GM2_global75={};
ActiveRxns_case1_GM1_local25={};
ActiveRxns_case1_GM2_local25={};
ActiveRxns_case1_GM1_local25_75={};
ActiveRxns_case1_GM2_local25_75={};
ActiveRxns_case1_GM1_local25_90={};
ActiveRxns_case1_GM2_local25_90={};
ActiveRxns_case2_GM1_global50={};
ActiveRxns_case2_GM2_global50={};
ActiveRxns_case2_GM1_global75={};
ActiveRxns_case2_GM2_global75={};
ActiveRxns_case2_GM1_local25={};
ActiveRxns_case2_GM2_local25={};
ActiveRxns_case2_GM1_local25_75={};
ActiveRxns_case2_GM2_local25_75={};
ActiveRxns_case2_GM1_local25_90={};
ActiveRxns_case2_GM2_local25_90={};

for i=1:length(expressionData.Tissue)
    ActiveRxns_case1_GM1_global50{i}=(find(expression.Rxns_case1_GM1(:,i)>=expression.ths_50));
    ActiveRxns_case1_GM2_global50{i}=(find(expression.Rxns_case1_GM2(:,i)>=expression.ths_50));
    ActiveRxns_case2_GM1_global50{i}=(find(expression.Rxns_case2_GM1_global50(:,i)>=expression.ths_50));
    ActiveRxns_case2_GM2_global50{i}=(find(expression.Rxns_case2_GM2_global50(:,i)>=expression.ths_50));
    ActiveRxns_case1_GM1_global75{i}=(find(expression.Rxns_case1_GM1(:,i)>=expression.ths_75));
    ActiveRxns_case1_GM2_global75{i}=(find(expression.Rxns_case1_GM2(:,i)>=expression.ths_75));
    ActiveRxns_case2_GM1_global75{i}=(find(expression.Rxns_case2_GM1_global75(:,i)>=expression.ths_75));
    ActiveRxns_case2_GM2_global75{i}=(find(expression.Rxns_case2_GM2_global75(:,i)>=expression.ths_75));
    ActiveRxns_case1_GM1_local25{i}=(find((expression.Rxns_case1_GM1(:,i)./expression.ths_case1_GM1_local25(:,i))>=1));
    ActiveRxns_case1_GM2_local25{i}=(find((expression.Rxns_case1_GM2(:,i)./expression.ths_case1_GM2_local25')>=1));
    ActiveRxns_case1_GM1_local25_75{i}=(find((expression.Rxns_case1_GM1(:,i)./expression.ths_case1_GM1_local25_75(:,i))>=1));
    ActiveRxns_case1_GM2_local25_75{i}=(find((expression.Rxns_case1_GM2(:,i)./expression.ths_case1_GM2_local25_75')>=1));
    ActiveRxns_case1_GM1_local25_90{i}=(find((expression.Rxns_case1_GM1(:,i)./expression.ths_case1_GM1_local25_90(:,i))>=1));
    ActiveRxns_case1_GM2_local25_90{i}=(find((expression.Rxns_case1_GM2(:,i)./expression.ths_case1_GM2_local25_90')>=1));
    ActiveRxns_case2_GM1_local25{i}=(find(expression.Rxns_case2_GM1_local25(:,i)>=1));
    ActiveRxns_case2_GM2_local25{i}=(find(expression.Rxns_case2_GM2_local25(:,i)>=1));
    ActiveRxns_case2_GM1_local25_75{i}=(find(expression.Rxns_case2_GM1_local25_75(:,i)>=1));
    ActiveRxns_case2_GM2_local25_75{i}=(find(expression.Rxns_case2_GM2_local25_75(:,i)>=1));
    ActiveRxns_case2_GM1_local25_90{i}=(find(expression.Rxns_case2_GM1_local25_90(:,i)>=1));
    ActiveRxns_case2_GM2_local25_90{i}=(find(expression.Rxns_case2_GM2_local25_90(:,i)>=1));
end

approaches=({'Case 1 GM1 Global T1 50th' ...
'Case 2 GM1 Global T1 50th' ...
'Case 1 GM2 Global T1 50th' ...
'Case 2 GM2 Global T1 50th' ...
'Case 1 GM1 Global T1 75th' ...
'Case 2 GM1 Global T1 75th' ...
'Case 1 GM2 Global T1 75th' ...
'Case 2 GM2 Global T1 75th' ...
'Case 1 GM1 Local T1 25th' ...
'Case 2 GM1 Local T1 25th' ...
'Case 1 GM2 Local T1 25th' ...
'Case 2 GM2 Local T1 25th' ...
'Case 1 GM1 Local T2 25th & 75th' ...
'Case 2 GM1 Local T2 25th & 75th' ...
'Case 1 GM2 Local T2 25th & 75th' ...
'Case 2 GM2 Local T2 25th & 75th' ...
'Case 1 GM1 Local T2 25th & 90th' ...
'Case 2 GM1 Local T2 25th & 90th' ...
'Case 1 GM2 Local T2 25th & 90th' ...
'Case 2 GM2 Local T2 25th & 90th'});

for k=1:20
    if k==1
        essGene=ActiveRxns_case1_GM1_global50;
    elseif k==2
        essGene=ActiveRxns_case2_GM1_global50;
    elseif k==3
        essGene=ActiveRxns_case1_GM2_global50;
    elseif k==4
        essGene=ActiveRxns_case2_GM2_global50;
    elseif k==5
        essGene=ActiveRxns_case1_GM1_global75;
    elseif k==6
        essGene=ActiveRxns_case2_GM1_global75;
    elseif k==7
        essGene=ActiveRxns_case1_GM2_global75;
   	elseif k==8
        essGene=ActiveRxns_case2_GM2_global75;
 	elseif k==9
        essGene=ActiveRxns_case1_GM1_local25;
   	elseif k==10
        essGene=ActiveRxns_case2_GM1_local25;
  	elseif k==11
        essGene=ActiveRxns_case1_GM2_local25;
 	elseif k==12
        essGene=ActiveRxns_case2_GM2_local25;
  	elseif k==13
        essGene=ActiveRxns_case1_GM1_local25_75;
 	elseif k==14
        essGene=ActiveRxns_case2_GM1_local25_75;
  	elseif k==15
        essGene=ActiveRxns_case1_GM2_local25_75;
 	elseif k==16
        essGene=ActiveRxns_case2_GM2_local25_75;
  	elseif k==17
        essGene=ActiveRxns_case1_GM1_local25_90;
 	elseif k==18
        essGene=ActiveRxns_case2_GM1_local25_90;
  	elseif k==19
        essGene=ActiveRxns_case1_GM2_local25_90;
  	elseif k==20
        essGene=ActiveRxns_case2_GM2_local25_90;
    end

    Jacc_index=zeros(length(essGene),length(essGene));
    for i=1:length(essGene)
        Gene1= (essGene{i});
        if isempty(Gene1)
            Jacc_index(i,:)= 0;
        else
            for j=1:length(essGene)
                Gene2=(essGene{j});
                if isempty(Gene1)
                    Jacc_index(i,j)=0;
                else
                    Jacc_index(i,j)= numel(intersect(Gene1,Gene2))/numel(union(Gene1,Gene2));
                end
            end
        end
    end
    figure(3);
        colormap hot
        subplot(5,4,k);imagesc(Jacc_index)
        title(approaches{k})
        cb = colorbar;
        ax=gca;
        ax.YTick=[1:length(expressionData.Tissue)];
        ax.YTickLabel={};
        ylabel('Tissues');
        ax.XTick=[1:length(expressionData.Tissue)];
        ax.XTickLabel={};
        xlabel('Tissues');
        caxis([0 1])
end

%% How many genes are considered as active depending on preprocessing combinations used
expression.countActiveGene_global_T1_25=[];
expression.countActiveGene_global_T1_50=[];
expression.countActiveGene_global_T1_75=[];
expression.countActiveGene_global_T1_90=[];
expression.countActiveGene_local_T1_25=[];
expression.countActiveGene_local_T2_25_75=[];
expression.countActiveGene_local_T2_25_90=[];
expression.ActiveGene_global_T1_25={};
expression.ActiveGene_global_T1_50={};
expression.ActiveGene_global_T1_75={};
expression.ActiveGene_global_T1_90={};
expression.ActiveGene_local_T1_25={};
expression.ActiveGene_local_T2_25_75={};
expression.ActiveGene_local_T2_25_90={};

for i=1:32
    expression.countActiveGene_global_T1_25(i)=length(find(expression.scoreGlobal(:,i)>=expression.ths_25));
    expression.countActiveGene_global_T1_50(i)=length(find(expression.scoreGlobal(:,i)>=expression.ths_50));
    expression.countActiveGene_global_T1_75(i)=length(find(expression.scoreGlobal(:,i)>=expression.ths_75));
    expression.countActiveGene_global_T1_90(i)=length(find(expression.scoreGlobal(:,i)>=expression.ths_90));
    expression.countActiveGene_local_T1_25(i)=length(find(expression.scoreLocal_T1_25(:,i)>=1));
    expression.countActiveGene_local_T2_25_75(i)=length(find(expression.scoreLocal_T2_25_75(:,i)>=1));
    expression.countActiveGene_local_T2_25_90(i)=length(find(expression.scoreLocal_T2_25_90(:,i)>=1));
    expression.ActiveGene_global_T1_25{i}=(find(expression.scoreGlobal(:,i)>=expression.ths_25));
    expression.ActiveGene_global_T1_50{i}=(find(expression.scoreGlobal(:,i)>=expression.ths_50));
    expression.ActiveGene_global_T1_75{i}=(find(expression.scoreGlobal(:,i)>=expression.ths_75));
    expression.ActiveGene_global_T1_90{i}=(find(expression.scoreGlobal(:,i)>=expression.ths_90));
    expression.ActiveGene_local_T1_25{i}=(find(expression.scoreLocal_T1_25(:,i)>=1));
    expression.ActiveGene_local_T2_25_75{i}=(find(expression.scoreLocal_T2_25_75(:,i)>=1));
    expression.ActiveGene_local_T2_25_90{i}=(find(expression.scoreLocal_T2_25_90(:,i)>=1));
end

up=[expression.countActiveGene_global_T1_25' ...
expression.countActiveGene_global_T1_50' ...
expression.countActiveGene_global_T1_75' ...
expression.countActiveGene_global_T1_90' ...
expression.countActiveGene_local_T1_25' ...
expression.countActiveGene_local_T2_25_75' ...
expression.countActiveGene_local_T2_25_90'];

up_percent=up./1663*100;
figure(4);boxplot(up_percent)
ax=gca;ax.XTick=1:7;
ax.XTickLabel=({'Global T1 25th','Global T1 50th','Global T1 75th','Global T1 90th','Local T1 25th','Local T2 25th & 75th','Local T2 25th & 90th'});
ax.XTickLabelRotation=(45);
ylabel('Percentage of genes considered as active','FontSize',14)

%% What is the similarity between genes considered as active depending on preprocessing combinations
approaches=({'Global T1 25th','Global T1 50th','Global T1 75th','Global T1 90th','Local T1 25th','Local T2 25th & 75th','Local T2 25th & 90th'});

for k=1:7
    if k==1
        essGene=expression.ActiveGene_global_T1_25;
    elseif k==2
        essGene=expression.ActiveGene_global_T1_50;
    elseif k==3
        essGene=expression.ActiveGene_global_T1_75;
    elseif k==4
        essGene=expression.ActiveGene_global_T1_90;
    elseif k==5
        essGene=expression.ActiveGene_local_T1_25;
    elseif k==6
        essGene=expression.ActiveGene_local_T2_25_75;
    elseif k==7
        essGene=expression.ActiveGene_local_T2_25_90;
    end
    Jacc_index=zeros(length(essGene),length(essGene));
    for i=1:length(essGene)
        Gene1= (essGene{i});
        if isempty(Gene1)
            Jacc_index(i,:)= 0;
        else
            for j=1:length(essGene)
                Gene2=(essGene{j});
                if isempty(Gene1)
                    Jacc_index(i,j)=0;
                else
                    Jacc_index(i,j)= numel(intersect(Gene1,Gene2))/numel(union(Gene1,Gene2));
                end
            end
        end
    end
    figure(5);
        colormap hot
        subplot(4,2,k);imagesc(Jacc_index)
        title(approaches{k})
        cb = colorbar;
        ax=gca;
        ax.YTick=([1:32]);ax.YTickLabel=({});ylabel('Tissues');
        ax.XTick=([1:32]);ax.XTickLabel=({});xlabel('Tissues');caxis([0 1])
end
