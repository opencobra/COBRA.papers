function  plotATPRxns(ATPcontribution, Allmodels,ind)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

new = {'ATPS4minew'
    'CYOOm2inew'
    'CYOOm3inew'
    'CYOR_u10minew'
    'NADH2_u10minew'
    'EX_cys_L(e)'};

original = {'ATPS4mi'
    'CYOOm2i'
    'CYOOm3i'
    'CYOR_u10mi'
    'NADH2_u10mi'};
if ~exist('ind','var')
    ind= 1; % default as ATP contribution
end
allmodelsValue=ATPcontribution;
% reorganize
if any(contains(fieldnames(allmodelsValue),{'SYNPD'}))
    names = [fieldnames(allmodelsValue.SYN); fieldnames(allmodelsValue.SYNPD)];
    allmodelsValue.allSYN = cell2struct([struct2cell(allmodelsValue.SYN); struct2cell(allmodelsValue.SYNPD)], names, 1);
    names = [fieldnames(allmodelsValue.ASYN); fieldnames(allmodelsValue.ASYNPD)];
    allmodelsValue.allASYN = cell2struct([struct2cell(allmodelsValue.ASYN); struct2cell(allmodelsValue.ASYNPD)], names, 1);
    allmodelsValue=rmfield(allmodelsValue,{'SYN','SYNPD','ASYN','ASYNPD'});
end
% for control and PD
groups=fieldnames(allmodelsValue);
for i=1:length(groups)
    subgroups=fieldnames(allmodelsValue.(groups{i}));
    if any(contains(subgroups,'constrain'))
        subgroups(contains(subgroups,'constrain'))=[];% ignore unconstrained models
    end
    if length(subgroups) == 4
        for j=1:2
            group=subgroups(contains(subgroups,string(j)));
            data1=allmodelsValue.(groups{i}).(group{1}).metRs;
            data1=sortrows(data1,2,'descend');
            data2=allmodelsValue.(groups{i}).(group{2}).metRs;
            data2=sortrows(data2,2,'descend');
            % identify all unique rxnIDs
            allrxns=unique([data1(:,1);data2(:,1)]);
            % re-order the rxns
            new1=cell(length(allrxns),3);
            for k=1:length(allrxns)
                if any(ismember(allrxns{k},data1(:,1)))
                    idx=find(ismember(data1(:,1),allrxns{k}));
                    new1(k,:)=data1(idx,:);
                else
                    new1(k,:)={allrxns{k},0,0};
                end
            end
            
            new1=sortrows(new1,2,'descend');
            allrxns=new1(:,1);
            % reorganize each data by allrxns
            new2=cell(length(allrxns),3);
            for k=1:length(allrxns)
                if any(ismember(allrxns{k},data2(:,1)))
                    idx=find(ismember(data2(:,1),allrxns{k}));
                    new2(k,:)=data2(idx,:);
                else
                    new2(k,:)={allrxns{k},0,0};
                end
            end
            
            d = figure('units','normalized','outerposition',[0 0 1 1]);
            fluxValue=[round(cell2mat(new1(:,2)),3),round(cell2mat(new2(:,2)),3)];% two group in one variable
            b1=bar(fluxValue);
            rxnsName=allrxns;
            % find subsystems
            %         index=ismember(models.(groups{i}).(subgroups{j}).rxns,rxnsName);
            %         subsystems=string(models.(groups{i}).(subgroups{j}).subSystems(index));
            %         xlabel= strcat(string(rxnsName), subsystems);
            %         xlabel=subsystems;
            xlabel=rxnsName;
            set(gca,'XTick',[1:length(xlabel)],'xticklabel',xlabel)
            set(gca,'XTickLabelRotation',25, 'FontSize',10)
            ylabel('Flux(umol/gDW/h)', 'FontSize',14,'FontWeight','bold')
            if isequal(ind,1)
                title(['Total ',char(string(length(allrxns))),' ATP contributing rxns in',' ',group{1}, ' and ',group{2}] ,'FontSize',12);
            else
                title(['Total ',char(string(length(allrxns))),' ATP consuming rxns in',' ',group{1}, ' and ',group{2}] ,'FontSize',12);
            end
            l1=char(string(length(data1(:,1))));
            l2=char(string(length(data2(:,1))));
            legend([group{1} '(total ' l1 ' ATP rxns)' ],[group{2} '(total ' l2 ' ATP rxns)' ]);
            legend('FontSize',14)
            legend('boxoff')
            
            for m=1:length(b1)
                text(b1(m).XEndPoints,b1(m).YEndPoints+double(m-0.9),string(fluxValue(:,m)), ...
                    'VerticalAlignment','bottom','horizontalalign','center')
            end
            
        end
    elseif length(subgroups) == 2 % only with SYN2 and SYNPD2
            data1=allmodelsValue.(groups{i}).(subgroups{1}).metRs;
            data1(cell2mat(data1(:,2))<1e-4,:)=[];
            data1=sortrows(data1,2,'descend');
            data2=allmodelsValue.(groups{i}).(subgroups{2}).metRs;
            data2(cell2mat(data2(:,2))<1e-4,:)=[];
            data2=sortrows(data2,2,'descend');
            % identify all unique rxnIDs
            allrxns=unique([data1(:,1);data2(:,1)]);
            % re-order the rxns
            new1=cell(length(allrxns),3);
            for k=1:length(allrxns)
                if any(ismember(allrxns{k},data1(:,1)))
                    idx=find(ismember(data1(:,1),allrxns{k}));
                    new1(k,:)=data1(idx,:);
                else
                    new1(k,:)={allrxns{k},0,0};
                end
            end
            
            new1=sortrows(new1,2,'descend');
            allrxns=new1(:,1);
            % reorganize each data by allrxns
            new2=cell(length(allrxns),3);
            for k=1:length(allrxns)
                if any(ismember(allrxns{k},data2(:,1)))
                    idx=find(ismember(data2(:,1),allrxns{k}));
                    new2(k,:)=data2(idx,:);
                else
                    new2(k,:)={allrxns{k},0,0};
                end
            end

            % Create a full-screen figure window
            %             d = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
            subplot(2,1,i)
            % Combine the flux values from two groups into one matrix and round the values
            fluxValue = round([cell2mat(new1(:,2)), cell2mat(new2(:,2))], 4);
            
            % Create a bar plot with a custom color scheme
            b1 = bar(fluxValue, 'FaceColor', 'flat');
            % Define a fashionable color palette
            colors = [0, 0.2941, 0.5843;  %  blue
                0.9255, 0.4784, 0.0314]; % orange
            
            for k = 1:size(fluxValue, 2)
                b1(k).CData = colors(k, :);
            end
            
            % Set the reaction names as x-axis labels
            rxnsName = allrxns;
            set(gca, 'XTick', 1:length(rxnsName), 'XTickLabel', rxnsName);
            
            % Customize x-axis label appearance
            set(gca, 'XTickLabelRotation', 25, 'FontSize', 10,'FontWeight', 'bold');
            
            % Remove the top x-axis ticks
            set(gca, 'Box', 'off');
            
            % Add a y-axis label
            ylabel('Flux (umol/gDW/h)', 'FontSize', 14, 'FontWeight', 'bold');
            ylim(gca,[0,max(max(fluxValue))+1])
            % Set the title based on whether reactions are contributing or consuming ATP
            if isequal(ind, 1)
                title(['Total ', num2str(length(allrxns)), ' ATP contribution reactions in ', subgroups{1}, ' and ', subgroups{2}], 'FontSize', 14);
            else
                title(['Total ', num2str(length(allrxns)), ' ATP consumption reactions in ', subgroups{1}, ' and ', subgroups{2}], 'FontSize', 14);
            end
            
            % Define the legend labels with the total number of ATP reactions
            l1 = num2str(length(data1(:,1)));
            l2 = num2str(length(data2(:,1)));

            if ind % 0=consumption; 1=contribution/production
                legend([subgroups{1}, ' (total ', l1, ' ATP generation reactions)'], [subgroups{2}, ' (total ', l2, ' ATP contribution reactions)'], ...
                    'FontSize', 12, 'Box', 'off');
            else
                legend([subgroups{1}, ' (total ', l1, ' ATP consumption reactions)'], [subgroups{2}, ' (total ', l2, ' ATP consumption reactions)'], ...
                    'FontSize', 12, 'Box', 'off');
            end
            
            % Annotate the bars with the flux values
            %             for m = 1:length(b1)
            %                 text(b1(m).XEndPoints, b1(m).YEndPoints + double(m - 0.9), string(fluxValue(:,m)), ...
            %                     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            %             end
            % Annotate the bars with the flux values
            text(b1(1).XEndPoints, b1(1).YEndPoints + double(0.1), string(fluxValue(:,1)), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            text(b1(2).XEndPoints, b1(2).YEndPoints + double(0.4), string(fluxValue(:,2)), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % add detailed rxns names
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            if any(contains(subgroups,'ASYN'))
                model=table();
                model.rxns=unique([Allmodels.ASYN.ASYN.rxns;Allmodels.ASYNPD.ASYNPD.rxns]);
                for k=1:length(model.rxns)
                    if ismember(model.rxns(k),Allmodels.ASYN.ASYN.rxns)
                        model.rxnNames(k)=Allmodels.ASYN.ASYN.rxnNames(ismember(Allmodels.ASYN.ASYN.rxns,model.rxns(k)));
                    else
                        model.rxnNames(k)=Allmodels.ASYNPD.ASYNPD.rxnNames(ismember(Allmodels.ASYNPD.ASYNPD.rxns,model.rxns(k)));
                    end
                end
            else
                model=table();
                model.rxns=unique([Allmodels.SYN.SYN.rxns;Allmodels.SYNPD.SYNPD.rxns]);
                for k=1:length(model.rxns)
                    if ismember(model.rxns(k),Allmodels.SYN.SYN.rxns)
                        model.rxnNames(k)=Allmodels.SYN.SYN.rxnNames(ismember(Allmodels.SYN.SYN.rxns,model.rxns(k)));
                    else
                        model.rxnNames(k)=Allmodels.SYNPD.SYNPD.rxnNames(ismember(Allmodels.SYNPD.SYNPD.rxns,model.rxns(k)));
                    end
                end
            end
            [A,B]=ismember(model.rxns,allrxns);
            B(B == 0) = [];
            Anno(B)=model.rxnNames(A);           
            Anno=Anno';
            annotationText= strcat(allrxns,' : ',Anno);
            
            % Number of columns and rows
            numCols = 2;
            numRows = ceil(length(annotationText) / numCols);
            
            % Define position settings
                columnWidth = 0.2; % Width of each column
                rowHeight = 0.01; % Height of each row
                marginX = 0.02; % Margin between columns
                marginY = 0.01; % Margin between rows
            % Create annotation boxes
            for idx = 1:length(annotationText)
                col = mod(idx-1, numCols); % Determine column
                row = floor((idx-1) / numCols); % Determine row
                
                % Calculate position
                xPos = 0.25 + col * (columnWidth + marginX);
                if i==1
                    yPos = 0.9 - row * (rowHeight + marginY);
                else
                    yPos = 0.42 - row * (rowHeight + marginY); % Adjusted for subplot
                end
                annotationPosition = [xPos, yPos, columnWidth, rowHeight];
                
                % Create the annotation box
                annotation('textbox', annotationPosition, 'String', annotationText{idx}, ...
                    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'none', ...
                    'FontSize', 8, 'HorizontalAlignment', 'left');
            end
            
            clear Anno annotationText
    end
    
    
end

