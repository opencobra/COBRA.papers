% CSF=readtable('~/drive/bioenergeticsPD/fromXi/data/csf_xomic.xlsx','sheet','mediaData',"TextType","string",'ReadVariableNames',true);
% dataFolder = ['~/drive/bioenergeticsPD/fromXi/reogranizeData'];
% bibliomicData = 'test_SYN.xlsx';
% specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep 'new3' filesep bibliomicData], 1, 1);
% specificData.mediaData(~ismember(specificData.mediaData.rxns,CSF.rxns),:)

% fad rxns in Recon3D
genericModelName = 'Recon3D_301.mat'; 
% genericModelName = 'Recon3DModel_301.mat'; % recon3DModel
switch genericModelName
    case 'Recon3DModel_301.mat'
        load('~/fork-cobratoolbox/test/models/mat/Recon3DModel_301.mat');
        modelOrig=model;
    case 'Recon3D_301.mat'
         load('~/fork-cobratoolbox/test/models/mat/Recon3D_301.mat');
        model=Recon3D;
end
%find fad in mito
% FADrxns=findRxnsFromMets(model,'fad[m]');
% formula=printRxnFormula(model,FADrxns);
% genes=model.grRules(ismember(model.rxns,FADrxns));
% FAD=[FADrxns,table(formula),genes];
% FAD.Properties.VariableNames={'FAD[m]rxns','formula','genes'};

if ismember('Recon3D_301.mat',genericModelName)
    modelOrig=model;
    FAD_changed=readtable('~/drive/bioenergeticsPD/fromXi/data/FADrxns/FAD_forRecon3D.xlsx','sheet','changed',"TextType","string",'ReadVariableNames',true);
    FAD_deleted=readtable('~/drive/bioenergeticsPD/fromXi/data/FADrxns/FAD_forRecon3D.xlsx','sheet','deleted',"TextType","string",'ReadVariableNames',true);
end

others=FAD(~ismember(FADrxns,FAD_changed.RxnID),:);
for i=1:length(others.("FAD[m]rxns"))
    others{i,3}=model.grRules(ismember(model.rxns,others.("FAD[m]rxns"){i}));
end
writetable(others,'~/Documents/others.xlsx')

genes2=model.grRules(ismember(model.rxns,FAD_changed.RxnID);
geneInfo=[table(FAD_changed.RxnID),table(genes2)];

writetable(geneInfo,'~/Documents/genes.xlsx')
%%
% all flavoprotein-related rxns in Recon3D
% load the human flavoprotein list from Agn's paper?
load('~/fork-cobratoolbox/test/models/mat/Recon3D_301.mat');
model=Recon3D;
flavoproteinAgn.Properties.VariableNames={'ProteinID','E_name','protein_name','Gene_name','EC_number','activity','GeneID','Disease','compartment','HGNC','cofactor'};
Genes=cellstr(num2cell(flavoproteinAgn.GeneID(:)));
flavoproteinAgn.Recon3D=ismember(Genes,model.genes);
rxns=findRxnsFromGenes(model,Genes);
genesRxn=fieldnames(rxns);
FlavoproteinRxns=[];
for i=1:length(fieldnames(rxns)) 
    rxn=rxns.(genesRxn{i});
    rxnID=rxn(:,1);
    FlavoproteinRxns=[FlavoproteinRxns;rxnID];
end
FlavoproteinRxns=unique(FlavoproteinRxns);
%
overlappedRxns=FlavoproteinRxns(ismember(FlavoproteinRxns,FAD_changed.RxnID));
unoverlappedRxns=FAD_changed.RxnID(~ismember(FAD_changed.RxnID,FADrxns));
unchangedFAD=FlavoproteinRxns(~ismember(FlavoproteinRxns,overlappedRxns));
%
% how many fad are free?
fads={'fad[c]','fad[m]','fad[n]','fad[g]','fad[r]','fad[x]','fad[i]','fad[l]','fad[e]'};
% allfadrxns=findRxnsFromMets(model,fads);
% [A,B]=find(model.S(ismember(model.mets,fads),:)<0); % fad redox rxns
% allfadrxns1=model.rxns(B);
[A,B]=find(model.S(ismember(model.mets,fads),:)); % all fad rxns
allfadrxns2=model.rxns(B);
fadrxnsWithoutFlavoprotein=allfadrxns2(~ismember(allfadrxns2,FlavoproteinRxns));
fadrxnsWithfadMets=allfadrxns2(ismember(allfadrxns2,FlavoproteinRxns));
FlavoproteinRxnsWithoutfadMets=FlavoproteinRxns(~ismember(FlavoproteinRxns,allfadrxns2));
%
% how many fmn are free?
fmns={'fmn[c]','fmn[m]','fmn[n]','fmn[g]','fmn[r]','fmn[x]','fmn[i]','fmn[l]','fmn[e]'};
[A,B]=find(model.S(ismember(model.mets,fmns),:)<0); % fad redox rxns
allfmnrxns1=model.rxns(B);
[A,B]=find(model.S(ismember(model.mets,fmns),:));
allfmnrxns2=model.rxns(B);
fmnrxnsWithoutFlavoprotein=allfmnrxns2(~ismember(allfmnrxns2,FlavoproteinRxns));
%% 
fads={'fad[c]','fad[m]','fad[n]','fad[g]','fad[r]','fad[x]','fad[i]','fad[l]','fad[e]'};
[A,B]=find(model.S(ismember(model.mets,fads),:)); % all fad rxns
allfadrxns2=model.rxns(B);
FAD_changed=readtable('~/drive/bioenergeticsPD/fromXi/data/FADrxns/FAD_forRecon3D.xlsx','sheet','changed',"TextType","string",'ReadVariableNames',true);
other = unique(allfadrxns2(~ismember(allfadrxns2,FAD_changed.RxnID)));
other = cell2table(other);
other.Properties.VariableNames={'RxnID'};
[A,B]=ismember(model.rxns,other.RxnID);
other.Subsystem(B(B~=0))=model.subSystems(A);
formula=printRxnFormula(model,other.RxnID);
other.Rawformula=formula;
other.GPR(B(B~=0))=model.grRules(A);
writetable(other,'~/Documents/other.xlsx')