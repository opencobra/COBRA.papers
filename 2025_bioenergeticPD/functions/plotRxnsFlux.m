function plotRxnsFlux(model1, model2,rxnList,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Running FBA
if ~isfield(param,'method')
    param.method='FBA';
end
switch param.method
    case 'FBA'
        changeCobraSolver('gurobi','all');
        FBAsolution1 = optimizeCbModel(model1,'max',1e-6);
        FBAsolution2 = optimizeCbModel(model2,'max',1e-6);
    case 'eFBA'
        param.solver = 'mosek';
        if strcmp(param.solver,'mosek')
            %set default mosek parameters for this type of problem
            param = mosekParamSetEFBA(param);
        end
        [FBAsolution1,~] = entropicFluxBalanceAnalysis(model1,param);
        [FBAsolution2,~] = entropicFluxBalanceAnalysis(model2,param);
end
%%
fluxValue=zeros(length(rxnList),2);
names=cell(length(rxnList),2);
names(:,1)=rxnList;
for i=1:length(rxnList)
    if find(ismember(rxnList{i},model1.rxns))
        index1(i)=find(ismember(model1.rxns,rxnList{i}));
        names{i,2}=model1.rxns(index1(i));
        fluxValue(i,1)=FBAsolution1.v(index1(i));
    else
        fluxValue(i,1)=0;
    end
    if find(ismember(rxnList{i},model2.rxns))
        index2(i)=find(ismember(model2.rxns,rxnList{i}));
        names{i,2}=model2.rxns(index2(i));
        fluxValue(i,2)=FBAsolution2.v(index2(i));
    else
        fluxValue(i,2)=0;
    end
    if isempty(names{i,2})
        names{i,2}=names{i,1};
    end
end

%%%%%%%%%%%%%       plot model1 and model2
            c = string(names(:,2));
            c = strrep(c,'_','-');
            m= size(c,1);
            figure
            b=bar (fluxValue)%,'FaceColor' , [0 0.4470 0.7410]
            set(gca,'XTick',[1:m],'xticklabel',c)
            set(gca,'XTickLabelRotation',20, 'FontSize',12)
            ylabel('Flux(umol/gDW)', 'FontSize',14,'FontWeight','bold')
            title(['rxn flux in two models with method ' param.method],'FontSize',12);
            legend('model1','model2');
            h = gca;
            h.YAxis.FontWeight = 'bold';
            h.XAxis.FontWeight = 'bold';
            hT=[];
            for i=1:length(b)
                text(b(i).XEndPoints,b(i).YEndPoints+0.5,string(fluxValue(:,i)), ...
                          'VerticalAlignment','bottom','horizontalalign','center')
            end
            legend('boxoff')
            
            
end

