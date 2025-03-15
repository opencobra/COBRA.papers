function plotATPproductionandConsumption(Allmodels,ATPproduction,ATPconsumption)
%UNTITLED Summary of this function goes here
%   plot the ATPproduction and ATPconsumption for each model 
group1=fieldnames(ATPproduction);
group2=fieldnames(ATPconsumption);

if length(group1)~=length(group2)
    error('The num of models in the ATPprodction must be same in the ATPconsumption')
end

for i=1:length(group1)
    subgroup1=fieldnames(ATPproduction.(group1{i}));
    subgroup1(contains(subgroup1,'constrain'))=[];% ignore unconstrained models
    for j=1:length(subgroup1)
        Production=ATPproduction.(group1{i}).(subgroup1{j}).metRs;
        Production=sortrows(Production,2,'descend');
        Consumption=ATPconsumption.(group1{i}).(subgroup1{j}).metRs;
        Consumption=sortrows(Consumption,2,'descend');
        d = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1)
        fluxValue=round(cell2mat(Production(:,2)),3);
        b1=bar(fluxValue);
        rxnsName=Production(:,1);
        for k=1:length(rxnsName)
            index{k}=find(ismember(Allmodels.(group1{i}).(subgroup1{j}).rxns,rxnsName{k}));
            subsystems{k}=string(Allmodels.(group1{i}).(subgroup1{j}).subSystems(index{k}));
        end
        subsystems=subsystems';
        %xlabel= strcat(string(rxnsName), subsystems);
        xlabel=subsystems;
        set(gca,'XTick',[1:length(xlabel)],'xticklabel',xlabel)
        set(gca,'XTickLabelRotation',25, 'FontSize',10)
        ylabel('Flux(umol/gDW/h)', 'FontSize',14,'FontWeight','bold')
        title(['ATP contributing/producing rxns in',' ',subgroup1{j}],'FontSize',12);
        x1=text(b1.XEndPoints,b1.YEndPoints+0.1,string(fluxValue), ...
                    'VerticalAlignment','bottom','horizontalalign','center');
        x2=text(b1.XEndPoints,b1.YEndPoints+1,string(rxnsName), ...
                    'VerticalAlignment','bottom','horizontalalign','center');
        set(x1, 'FontSize',8)
        set(x2, 'FontSize',8)
        %set(x2,'Rotation',15, 'FontSize',8)  
            
        subplot(1,2,2)
        b2=bar(cell2mat(Consumption(:,2)));
        fluxValue=round(cell2mat(Consumption(:,2)),3);
        rxnsName=Consumption(:,1);
        clear index subsystems
        for k=1:length(rxnsName)
            index{k}=find(ismember(Allmodels.(group1{i}).(subgroup1{j}).rxns,rxnsName{k}));
            subsystems{k}=string(Allmodels.(group1{i}).(subgroup1{j}).subSystems(index{k}));
        end
        subsystems=subsystems';
        %xlabel= strcat(string(rxnsName), subsystems);
        xlabel=subsystems;
        set(gca,'XTick',[1:length(xlabel)],'xticklabel',xlabel)
        set(gca,'XTickLabelRotation',25, 'FontSize',10)
        ylabel('Flux(umol/gDW/h)', 'FontSize',14,'FontWeight','bold')
        title(['ATP consumption rxns in',' ',subgroup1{j}],'FontSize',12);
        x1=text(b2.XEndPoints,b2.YEndPoints+0.3,string(fluxValue), ...
                    'VerticalAlignment','bottom','horizontalalign','center');
        set(x1, 'FontSize',8)
        x2=text(b2.XEndPoints,b2.YEndPoints+2,string(rxnsName), ...
                    'VerticalAlignment','bottom','horizontalalign','center');
        set(x2, 'FontSize',8)
    end
end

end

