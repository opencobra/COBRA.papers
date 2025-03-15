%% input data folder
inputFolder = '~/drive/metaPD/data/metabolomics/inputData';
% specificData
%bibliomicData = 'PD_alldata.xlsx'; %replicated mets + non replicated mets
bibliomicData = 'PD_test4_withcorerxns.xlsx';%all replicated mets+corerxns

%% Select the generic model name
genericModelName = 'modelDecompSCT_new2';

%% Select the results folder
resultsFolder = ['~/drive/metaPD/results/metabolomics/model/' genericModelName];
if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder)
end

%% Load generic model
switch genericModelName
    case 'recon3Dmodel'
        load('Recon3DModel_301.mat');
    case 'reconX'
        load('ReconX.mat')
    case 'modelDecomp'
        %load a decompartmentalised Recon
        load('~/drive/metaPD/results/metabolomics/model/modelDecomp.mat')
    case 'modelDecomp_SC'
        %% load  a decompartmentalised SConsistentsubset Recon
        load('~/drive/metaPD/results/metabolomics/model/modelDecomp_SC.mat')
        model = rmfield(model,{'SInConsistentMetBool','SInConsistentRxnBool','SConsistentMetBool','SConsistentRxnBool'});%'fluxConsistentMetBool','fluxConsistentRxnBool','fluxInConsistentMetBool','fluxInConsistentRxnBool'
    case 'modelDecompSCT'
        %% load  a decompartmentalised stoichiometrically, flux and thermodynamically flux consistent subset of Recon
        load('~/drive/metaPD/results/metabolomics/model/20221215T094823_modelDecompSCT.mat')
    case 'modelDecompSCT_new1'
        % load  a decompartmentalised stoichiometrically, flux and
        % thermodynamically flux consistent subset of Recon(only remove Atifical rxns)
        load('~/drive/metaPD/results/metabolomics/model/modelDecompSCT_new1.mat')
    case 'modelDecompSCT_new2'
        % load  a decompartmentalised stoichiometrically, flux and
        % thermodynamically flux consistent subset of Recon(remove Atifical rxns ,'R group synthesis','Pool reactions')
        load('~/drive/metaPD/results/metabolomics/model/modelDecompSCT_new2.mat')
    case 'modelDecompSCT_new3'
        % load  a decompartmentalised stoichiometrically, flux and
        % thermodynamically flux consistent subset of Recon(remove Atifical rxns,'R group synthesis','Pool reactions','Protein formation','Protein modification','Protein assembly')
        load('~/drive/metaPD/results/metabolomics/model/modelDecompSCT_new3.mat')
end

if 1
    %% Identify highly connected metabolites
%     param.n = 1000; % Connectivity of top x metabolites
    param.plot = 0; % Do not plot ranked connectivity
    param.internal = 1; % Ignore connectivity of stoichiometrically inconsistent part
    [rankMetConnectivity,rankMetInd,rankConnectivity] = rankMetabolicConnectivity(model, param);
    param.n = sum(rankConnectivity>9);
    boolConnected = false(size(model.S,1),1);
    boolConnected(rankMetInd(1:param.n))=1;
    
    %set weights corresponding to the highest connected metabolites to NaN so they are relaxed, unless they are core metabolites
    metWeights = zeros(size(model.S,1),1);
    metWeights(boolConnected)=NaN;
else
    metWeights = zeros(size(model.S,1),1);
end


%% Prepare data for xomics
specificData = preprocessingOmicsModel([inputFolder filesep bibliomicData], 1, 1);

%these weights were too small I think should be minimum -1 for a metabolite or reaction that must be present
specificData.presentMetabolites.weights = specificData.presentMetabolites.weights*100;

%% negative weights from clinical metabolomic data
metWeights = mapAontoB(specificData.presentMetabolites.mets,model.mets,specificData.presentMetabolites.weights,metWeights);%-specificData.presentMetabolites.Entropy

%% small penalty on other metabolites
metWeights(metWeights==0)=0.01;

%% small penaly on reaction weights
rxnWeights = 0.01*ones(size(model.S,2),1);

%% Extract a model that satisfies steady state constraints only where metWeights is not NaN
for i=1:10
param.plotThermoKernelStats=0;
param.plotThermoKernelWeights=0;     
param.printLevel = 1;
param.findThermoConsistentFluxSubset = 0;
param.saveModelSFC = 0;
[thermoModel, thermoModelMetBool, thermoModelRxnBool] = extractSemiSteadyStateModel(model,rxnWeights, metWeights, param);
name=strcat(['model', char(string((i)))]);
multimodels.(name)=thermoModel;
end
%% 
[overlapResults,statistic]=compareXomicsModels(multimodels);
plotOverlapResults(overlapResults,statistic);
%%
corerxn=ismember(model.rxns,overlapResults.rxns.alloverlap);
rxnWeights = 0.01*ones(size(model.S,2),1);
rxnWeights(corerxn)=-1.1;
%% find a stable rxnWeights from overlapped part
%contain all the overlapped rxns in rxnWeight (add corerxns1 over the corerxns)
n=0;
go=1;
while go==1
    rxnWeights1=rxnWeights;
    for j=1:10
        param.plotThermoKernelStats=0;
        param.plotThermoKernelWeights=0;
        param.printLevel = 1;
        param.findThermoConsistentFluxSubset = 0;
        param.saveModelSFC = 0;
        [thermoModel, thermoModelMetBool, thermoModelRxnBool] = extractSemiSteadyStateModel(model,rxnWeights, metWeights, param);
        name=strcat(['model', char(string((j)))]);
        multimodels.(name)=thermoModel;
    end
    
    [overlapResults,statistic]=compareXomicsModels(multimodels);
    corerxn=ismember(model.rxns,overlapResults.rxns.alloverlap);
    rxnWeights(corerxn)=-1.1;
    n=n+1
    
    percentage=length(overlapResults.rxns.alloverlap)./max(statistic.overlapnumber_rxns{:,2:end});
    
    if ~any(rxnWeights1~=rxnWeights) & percentage>=0.85
        go=0;
    end
end
%% get the overlapped plot
[overlapResults,statistic]=compareXomicsModels(multimodels);
plotOverlapResults(overlapResults,statistic)
%%
allcorerxns=findRxnsFromMets(model,specificData.presentMetabolites.mets);
coremodel=extractSubNetwork(model,model.rxns(rxnWeights==-1.1));
coremodel = updateGenes(coremodel);
corerxns=allcorerxns(ismember(allcorerxns,rxns));
corecoremodel=extractSubNetwork(coremodel,corerxns);
% subsystems
subs=tabulate(string(coremodel.subSystems));
subs=sortrows(subs,2,'d');
%change compartment and rxnsName to draw the Ecsher map
coremodel.rxnNames=coremodel.rxns;
coremodel.mets=strcat(coremodel.mets,'[no]');
corecoremodel.rxnNames=corecoremodel.rxns;
corecoremodel.mets=strcat(corecoremodel.mets,'[no]');
%% metsData for Escher map
load('~/drive/metaPD/results/metabolomics/data/v2/diagnosis/subset/metaPD.mat')
data=metaPD.total(ismember(metaPD.total.all,coremodel.mets),:);
for i=1:(length(data.all))
    if data.inconsistent(i)
        data.metsdata(i)=3;
    end
    if data.inconsistent(i) & (data.increased_highfreq(i) | data.decreased_highfreq(i))
        data.metsdata(i)=4;
    end
    if data.inconsistent(i)==0 & data.increased_Realfrequency(i)==1;
        data.metsdata(i)=2;
    end
     if data.inconsistent(i)==0 & data.decreased_Realfrequency(i)==1;
        data.metsdata(i)=1;
    end
     if data.inconsistent(i)==0 & data.increased_Realfrequency(i)>1;
        data.metsdata(i)=6;
    end
     if data.inconsistent(i)==0 & data.decreased_Realfrequency(i)>1;
        data.metsdata(i)=5;
    end   
end
%add caffeine
data=data(:,[1,size(data,2)]);
cafe=specificData.presentMetabolites([10:12,21,25,38:39,end-2:end],:);
cafemetsdata=[5;5;5;4;4;5;5;5;5;5];
cafemetsdata=table(cafe.mets,cafemetsdata);
data.Properties.VariableNames= {'metabolties';'value'};
cafemetsdata.Properties.VariableNames= {'metabolties';'value'};
%
data=[data ; cafemetsdata];
%data.metabolties=strcat(data.metabolties,'[no]');
%writetable(data,'~/drive/metaPD/results/metabolomics/model/modelDecompSCT_new2/164/coremodel/formap/mapfile/metsData.csv')
%% top subs
%load('modelDecompSCT_new2.mat')
allsubs=tabulate(string(model.subSystems));
allsubs=sortrows(allsubs,2,'d');
allsubs=cell2table(allsubs);
% subsystems
subs=tabulate(string(coremodel.subSystems));
subs=sortrows(subs,2,'d');
subs=cell2table(subs);
subs.subs3=[];
for i=1:length(subs.subs1)
subs.subs3(i)=allsubs.allsubs2(ismember(allsubs.allsubs1,subs.subs1(i)));
end
subs.RxnRatio=round(subs.subs2./subs.subs3,2);

for i=1:length(subs.subs1)
    rxns=findRxnsFromSubSystem(coremodel,subs.subs1(i));
    mets=findMetsFromRxns(coremodel,rxns);
    subs.CoreMet(i)=sum(ismember(mets,specificData.presentMetabolites.mets));
    allrxns=findRxnsFromSubSystem(model,subs.subs1(i));
    allmets=findMetsFromRxns(model,allrxns);
    subs.AllCoreMets(i)=sum(ismember(allmets,specificData.presentMetabolites.mets));
    subs.CoreMetRatio(i)=round(subs.CoreMet(i)/subs.AllCoreMets(i),4);
    subs.allmets(i)=sum(ismember(mets,data.metabolties));
    subs.AllMetRatio(i)=round(subs.allmets(i)./217,4);
end

subs.Properties.VariableNames={'Subsystems','rxnNum','AllRxnNum','RxnRatio','CoreMets','AllCoreMets','CoreMetRatio','AllMets','AllMetRatio'};

ignore={'Exchange/demand reaction','Protein formation','Miscellaneous','Transport, extracellular','Peptide metabolism'};

subs(ismember(subs.Subsystems,ignore),:)=[];
subs=sortrows(subs,5,'d');
%% get top 20 subsystems 
topsubs=subs(1:20,:);
topsubs.CoreMetRatio=round(topsubs.CoreMetRatio,2);
topsubs=sortrows(topsubs,5,'d'); %ordered subsystems by the number of core metabolites
%% unique core mets fraction
for i=1:length(topsubs.Subsystems)
    rxns=findRxnsFromSubSystem(coremodel,topsubs.Subsystems(i));
    mets=findMetsFromRxns(coremodel,rxns);
    subname=subs.Subsystems{i};
    CoreMet{i,1}=mets(ismember(mets,specificData.presentMetabolites.mets));
end

for i=1:20
    CoreMet{i,2}=CoreMet{i,1};
    for j=1:20
        if i~=j
            CoreMet{i,2}=CoreMet{i,2}(~ismember(CoreMet{i,2},CoreMet{j,1}));
        else
            continue
        end
    end
end
for i=1:20
CoreMet{i,3}=size(CoreMet{i,2},1)/topsubs.CoreMets(i); 
end
%% exclusive mets in each subsystem
for i=1:length(subs.Subsystems)
    colBool=ismember(coremodel.rxns,findRxnsFromSubSystem(coremodel,subs.Subsystems(i)));
    rowBool=ones(length(coremodel.mets),1);
    colBool=logical(colBool);
    rowBool=logical(rowBool);
    restricedRowBool = getCorrespondingRows(coremodel.S, rowBool, colBool, 'exclusive');
    subs.AlluniqueMets(i)=size(coremodel.mets(restricedRowBool),1);
    subs.CoreUniqueMets(i)=sum(ismember(specificData.presentMetabolites.mets,coremodel.mets(restricedRowBool)));
    subs.CoreUniqueMetsRatio(i)=subs.CoreUniqueMets(i)/subs.AlluniqueMets(i);
end
%% subsystem statistics
randomData=array2table(zeros(20,100));
randomData.Properties.RowNames=topsubs.Subsystems;
for i=1:100
    % Get the total number of reactions in the model
    numReactions = length(model.rxns);
    % Set the number of reactions you want to randomly pick
    numRandomReactions = length(coremodel.rxns);
    
    % Generate a random index vector of size numRandomReactions
    randomIndex = randperm(numReactions, numRandomReactions);
    
    % Get the reactions corresponding to the random indices
    randomReactions = model.rxns(randomIndex);
    
    % generate the submodel
    name=strcat(['model', char(string((i)))]);
    randomModels.(name)=extractSubNetwork(model,randomReactions);
    for j=1:length(topsubs.Subsystems)
        randomData{j,i}=size(findRxnsFromSubSystem(randomModels.(name),topsubs.Subsystems(j)),1);
    end
end

for i=1:length(topsubs.Subsystems)
    % Input data
    data = randomData(i,:);
    data=table2array(data);
    x = table2array(topsubs(i,2));% the value for which you want to calculate the p-value
    
    % Calculate t statistic and degrees of freedom
    n = size(data,2);
    df = n - 1;
    mean_data = mean(data);
    topsubs.mean(i)=mean_data;
    std_data = std(data);
    topsubs.std(i)=std_data;
    
    t_stat = (mean_data - x) / (std_data / sqrt(n));
    % Calculate 2-tailed p-value using the cumulative distribution function (CDF) of the t-distribution
    p_value = 2*(1-tcdf(abs(t_stat),df));
    topsubs.PValue(i)=p_value;
    
    %disp(['The p-value for x = ', num2str(x), ' is ', num2str(p_value)])
end
%% add bar chart and box plot for subsystems
%
figure('Units', 'pixels', 'Position', [100, 100, 1000, 700])
hold on
topsubs.CoreMetProp=topsubs.CoreMets./137;
b1=barh(-1*topsubs.RxnRatio);
b2=barh(topsubs.CoreMetProp);
for i = 1:20
    boxplot(-1*randomData{i,:}/topsubs.AllRxnNum(i),i,'orientation','horizontal','positions',i, 'Whisker', 1,'Symbol','');
end
hold off
%set(gca, 'box', 'off')
labels=string(topsubs.Subsystems);
set(gca,'Ydir','reverse')
yti=1:20;
set(gca,'ytick',yti,'yticklabel',labels(1:1:end));
yticks(1:20);
%set(gca,'ytick',yti);
%b.FaceColor="#D95319";
%add bar labels
%for rxns
xtips1=b1.YEndPoints;
ytips1=b1.XEndPoints;
%labels=string(topsubs.rxnNum);
%text(xtips1-0.01, ytips1,labels,'VerticalAlignment','middle','Horiz','right');
%for mets
xtips2=b2.YEndPoints;
ytips2=b2.XEndPoints;
%labels=string(topsubs.CoreMets);
%text(xtips2-0.1, ytips2,labels,'VerticalAlignment','middle','Horiz','left');
% for rxns
labels=string(topsubs.RxnRatio);
text(xtips1-0.001, ytips1,labels,'VerticalAlignment','middle','Horiz','right');
% for mets
labels=string(round(topsubs.CoreMetProp,2));
text(xtips2+0.005, ytips2,labels,'VerticalAlignment','middle','Horiz','left');

xlabel('Reaction fraction vs CoreMetabolites fraction')
ylabel('Subsystems')
title('Top 20 subsystems in the model','VerticalAlignment','bottom')
xticks([])
ax=gca;
ax.YAxisLocation='left';
ax.TickDir='out'
xlim([-1.2,0.26])

% Add *** at end of chart
hold on;
text(xtips1-0.23, ytips1, '*','VerticalAlignment','middle', 'HorizontalAlignment', 'right', 'FontSize', 14);
hold off;

% Create the legend
legend( {'Reactions', 'CoreMetabolites'},'Location','best');
legend('boxoff')


