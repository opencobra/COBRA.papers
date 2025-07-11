% Calculate ATP contribution with different ATPM
types=fieldnames(Allmodels);
for i=1:length(types)
    %define the samples:
    models = fieldnames(Allmodels.(types{i}));
    models=models(~contains(models,'constrain'));% ignore unconstrained models
    for j=1:length(models)
        switch range
            case '10-100'
                if contains(models{j},'ASTRO')
                    d = 1:10:50; %reduced to the maximum ATPM in astrocytes
                else
                    %             d=10:50:700; % for both SYN and ASYN
                    d=10:10:100;
                end
            case '10-660'
                if contains(models{j},'ASTRO')
                    d = 1:10:50; %reduced to the maximum ATPM in astrocytes
                else
                    d=10:50:700;
                end
        end     
        demandATP.(types{i}).(models{j}).atpm = zeros(length(d),1);
        demandATP.(types{i}).(models{j}).rxns = {};
        demandATP.(types{i}).(models{j}).total = {};
        demandATP.(types{i}).(models{j}).highest = {};
        demandATP.(types{i}).(models{j}).glycolysis = {};
        for k = 1:length(d)
            demand = d(1,k);
            model_FVA = Allmodels.(types{i}).(models{j});
            model_FVA = changeRxnBounds(model_FVA, 'ATPM', demand, 'b');
%             method='FBA';
            method='eFBA';
            eucNorm = 1e-4;
            genericMetaboContribution;
            %flux proportion
            demandATP.(types{i}).(models{j}).atpm(k,1) = d(1,k);
            demandATP.(types{i}).(models{j}).rxns{k,1} = metaboContribution.metabo.reactions;%all ATP production rxns
            demandATP.(types{i}).(models{j}).total{k,1} = sum([metaboContribution.metabo.reactions{:,2}]);%all atp production flux
            demandATP.(types{i}).(models{j}).highest{k,1} = metaboContribution.metabo.maxContRxn{1,1};%highest contribution rxn
            demandATP.(types{i}).(models{j}).PGK(k,1) = metaboPGK; %PGK flux proportion
            demandATP.(types{i}).(models{j}).PYK(k,1) = metaboPYK;
            demandATP.(types{i}).(models{j}).ATPS4mi(k,1) = metaboATPS4mi;
            demandATP.(types{i}).(models{j}).r0408(k,1) = metabor0408;
            demandATP.(types{i}).(models{j}).SUCOASm(k,1) = metaboSUCOASm;
            demandATP.(types{i}).(models{j}).NDPK6(k,1) = metaboNDPK6;
            demandATP.(types{i}).(models{j}).RE0124C(k,1) = metaboRE0124C;
            demandATP.(types{i}).(models{j}).NDPK2(k,1) = metaboNDPK2;
            demandATP.(types{i}).(models{j}).NMNATr(k,1) = metaboNMNATr;
            demandATP.(types{i}).(models{j}).glycolysis{k,1} = cell2mat(metaboPGK) + cell2mat(metaboPYK); %PGK+PYK proportion
            %flux value
            demandATP.(types{i}).(models{j}).PGK(k,2) = metaboPGKvalue; 
            demandATP.(types{i}).(models{j}).PYK(k,2) = metaboPYKvalue;
            demandATP.(types{i}).(models{j}).ATPS4mi(k,2) = metaboATPS4mivalue;
            demandATP.(types{i}).(models{j}).r0408(k,2) = metabor0408value;
            demandATP.(types{i}).(models{j}).SUCOASm(k,2) = metaboSUCOASmvalue;
            demandATP.(types{i}).(models{j}).NDPK6(k,2) = metaboNDPK6value;
            demandATP.(types{i}).(models{j}).RE0124C(k,2) = metaboRE0124Cvalue;
            demandATP.(types{i}).(models{j}).NDPK2(k,2) = metaboNDPK2value;
            demandATP.(types{i}).(models{j}).NMNATr(k,2) = metaboNMNATrvalue;
            demandATP.(types{i}).(models{j}).glycolysis{k,2} = cell2mat(metaboPGKvalue) + cell2mat(metaboPYKvalue); %PGK+PYK flux value
        end
    end
end
%%
% add flux proportion for UMPK and URIDK3;
% ADK1;NDPK1;NDPK3;NDPK4;NDPK5;NDPK8;NDPK9;NDPK10;r0345
types=fieldnames(demandATP);
for i=1:length(types)
    models=fieldnames(demandATP.(types{i}));
    for j=1:length(models)
        modelATPrxns=demandATP.(types{i}).(models{j});
        for k=1:length(modelATPrxns.rxns)
            ATPrxns=modelATPrxns.rxns{k};
            if any(ismember('UMPK',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).UMPK{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'UMPK'),3))*100;
                demandATP.(types{i}).(models{j}).UMPK{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'UMPK'),2));
            else
                demandATP.(types{i}).(models{j}).UMPK{k,1}=0;
                demandATP.(types{i}).(models{j}).UMPK{k,2}=0;
            end
            if any(ismember('URIDK3',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).URIDK3{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'URIDK3'),3))*100;
                demandATP.(types{i}).(models{j}).URIDK3{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'URIDK3'),2));
            else
                demandATP.(types{i}).(models{j}).URIDK3{k,1}=0;
                demandATP.(types{i}).(models{j}).URIDK3{k,2}=0;
            end
            if any(ismember('ADK1',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).ADK1{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'ADK1'),3))*100;
                demandATP.(types{i}).(models{j}).ADK1{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'ADK1'),2));
            else
                demandATP.(types{i}).(models{j}).ADK1{k,1}=0;
                demandATP.(types{i}).(models{j}).ADK1{k,2}=0;
            end
            if any(ismember('NDPK1',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).NDPK1{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK1'),3))*100;
                demandATP.(types{i}).(models{j}).NDPK1{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK1'),2));
            else
                demandATP.(types{i}).(models{j}).NDPK1{k,1}=0;
                demandATP.(types{i}).(models{j}).NDPK1{k,2}=0;
            end
            if any(ismember('NDPK3',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).NDPK3{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK3'),3))*100;
                demandATP.(types{i}).(models{j}).NDPK3{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK3'),2));
            else
                demandATP.(types{i}).(models{j}).NDPK3{k,1}=0;
                demandATP.(types{i}).(models{j}).NDPK3{k,2}=0;
            end
            if any(ismember('NDPK4',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).NDPK4{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK4'),3))*100;
                demandATP.(types{i}).(models{j}).NDPK4{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK4'),2));
            else
                demandATP.(types{i}).(models{j}).NDPK4{k,1}=0;
                demandATP.(types{i}).(models{j}).NDPK4{k,2}=0;
            end
            if any(ismember('NDPK5',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).NDPK5{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK5'),3))*100;
                demandATP.(types{i}).(models{j}).NDPK5{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK5'),2));
            else
                demandATP.(types{i}).(models{j}).NDPK5{k,1}=0;
                demandATP.(types{i}).(models{j}).NDPK5{k,2}=0;
            end
            if any(ismember('NDPK7',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).NDPK7{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK7'),3))*100;
                demandATP.(types{i}).(models{j}).NDPK7{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK7'),2));
            else
                demandATP.(types{i}).(models{j}).NDPK7{k,1}=0;
                demandATP.(types{i}).(models{j}).NDPK7{k,2}=0;
            end
            if any(ismember('NDPK8',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).NDPK8{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK8'),3))*100;
                demandATP.(types{i}).(models{j}).NDPK8{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK8'),2));
            else
                demandATP.(types{i}).(models{j}).NDPK8{k,1}=0;
                demandATP.(types{i}).(models{j}).NDPK8{k,2}=0;
            end
            if any(ismember('NDPK9',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).NDPK9{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK9'),3))*100;
                demandATP.(types{i}).(models{j}).NDPK9{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK9'),2));
            else
                demandATP.(types{i}).(models{j}).NDPK9{k,1}=0;
                demandATP.(types{i}).(models{j}).NDPK9{k,2}=0;
            end
            if any(ismember('NDPK10',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).NDPK10{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK10'),3))*100;
                demandATP.(types{i}).(models{j}).NDPK10{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'NDPK10'),2));
            else
                demandATP.(types{i}).(models{j}).NDPK10{k,1}=0;
                demandATP.(types{i}).(models{j}).NDPK10{k,2}=0;
            end
            if any(ismember('r0345',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).r0345{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'r0345'),3))*100;
                 demandATP.(types{i}).(models{j}).r0345{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'r0345'),2));
            else
                demandATP.(types{i}).(models{j}).r0345{k,1}=0;
                demandATP.(types{i}).(models{j}).r0345{k,2}=0;
            end
            if any(ismember('CYTK1',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).CYTK1{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'UMPK'),3))*100;
                demandATP.(types{i}).(models{j}).CYTK1{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'UMPK'),2));
            else
                demandATP.(types{i}).(models{j}).CYTK1{k,1}=0;
                demandATP.(types{i}).(models{j}).CYTK1{k,2}=0;
            end
            if any(ismember('CYTK2',ATPrxns(:,1)))
                demandATP.(types{i}).(models{j}).CYTK2{k,1}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'UMPK'),3))*100;
                demandATP.(types{i}).(models{j}).CYTK2{k,2}=cell2mat(ATPrxns(ismember(ATPrxns(:,1),'UMPK'),2));
            else
                demandATP.(types{i}).(models{j}).CYTK2{k,1}=0;
                demandATP.(types{i}).(models{j}).CYTK2{k,2}=0;
            end
        end
    end
end