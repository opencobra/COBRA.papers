%perform Sanity Checks
types=fieldnames(Allmodels);
% method='FBA';
method='eFBA';
runSingleGeneDeletion=0;
extraCellCompIn = '[e]';
extraCellCompOut = '[e]';
for i=1:length(types)
    %define the models:
    models = fieldnames(Allmodels.(types{i}));
    models=models(~contains(models,'constrain'));% ignore unconstrained models
    for j=1:length(models)
        model=Allmodels.(types{i}).(models{j});
        [TableChecks.(types{i}).(models{j}), Table_csources.(types{i}).(models{j}), CSourcesTestedRxns.(types{i}).(models{j}), TestSolutionNameOpenSinks.(types{i}).(models{j}),~] = performSanityChecks_bioenergeticsPD(model,...
            extraCellCompIn,extraCellCompOut,runSingleGeneDeletion,method);
    end
end
%%
types=fieldnames(TestSolutionNameOpenSinks);
k=1;
for i=1:length(types)
    models=fieldnames(TestSolutionNameOpenSinks.(types{i}));
    for j=1:length(models)
        eachmodel_testHumanFct=TestSolutionNameOpenSinks.(types{i}).(models{j});
        eachtable=TableChecks.(types{i}).(models{j});
%         eachsource=Table_csources.(types{i}).(models{j});
        %eachmodel_tablechecks=TableChecks.(types{i}).(models{j});
        AllTestSolution(:,k)=eachmodel_testHumanFct(:,2);
        TableCheck(:,k)=eachtable(:,2);
%         TableCsources=[TableCsources,eachsource(:,[2:3])];
        k=k+1;
    end
end

all=size(AllTestSolution,2);
for m=1:size(AllTestSolution,1)
    eachrow=string(AllTestSolution(m,1:2));
    if any(eachrow=='0')
        index=eachrow=='0';
        eachrow(index)=NaN;
    end
    
    if any(~ismissing(eachrow))
        AllTestSolution{m,all+1} = 'YES';
    else
        AllTestSolution{m,all+1} = 'NO';
    end
end

index=string(AllTestSolution(:,end))=='NO';

AllTestSolution=cell2table(AllTestSolution);
% if isempty(eachmodel_testHumanFct{324,1})
%     eachmodel_testHumanFct(324,1)={'uptake of cholic acid - CHOLATEt3'};
% end
AllTestSolution.Properties.RowNames=eachmodel_testHumanFct(:,1);
if all==20
    AllTestSolution.Properties.VariableNames={'SYN1','SYN1Uncon','SYN2','SYN2Uncon','SYNPD1','SYNPD1Uncon','SYNPD2','SYNPD2Uncon',...
        'ASYN1','ASYN1Uncon','ASYN2','ASYN2Uncon','ASYNPD1','ASYNPD1Uncon','ASYNPD2','ASYNPD2Uncon','Astro1','Astro1Uncon','Astro2','Astro2Uncon','index'};
elseif all==10
    AllTestSolution.Properties.VariableNames={'SYN1','SYN2','SYNPD1','SYNPD2','ASYN1','ASYN2','ASYNPD1','ASYNPD2','Astro1','Astro2','index'};
elseif all == 4
    AllTestSolution.Properties.VariableNames={'SYN2','SYNPD2','ASYN2','ASYNPD2','index'};
end

AllTestSolution(index,:)=[];
%
TableCheck=[eachtable(:,1),TableCheck];