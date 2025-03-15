function plots_CandPD_v2(demandATP,f,d)
%UNTITLED3 Summary of this function goes here
%   change demandATP structure
if ~exist('f','var')
    f='barchart';
end

if ~exist('d','var')
    error('check input data');
end
if any(contains(fieldnames(demandATP),{'SYNPD'})) & any(contains(fieldnames(demandATP),{'ASYNPD'}))
    names = [fieldnames(demandATP.SYN); fieldnames(demandATP.SYNPD)];
    demandATP.allSYN = cell2struct([struct2cell(demandATP.SYN); struct2cell(demandATP.SYNPD)], names, 1);
    names = [fieldnames(demandATP.ASYN); fieldnames(demandATP.ASYNPD)];
    demandATP.allASYN = cell2struct([struct2cell(demandATP.ASYN); struct2cell(demandATP.ASYNPD)], names, 1);
    demandATP=rmfield(demandATP,{'SYN','SYNPD','ASYN','ASYNPD'});
elseif any(contains(fieldnames(demandATP),{'SYNPD'}))
    names = [fieldnames(demandATP.SYN); fieldnames(demandATP.SYNPD)];
    demandATP.allSYN = cell2struct([struct2cell(demandATP.SYN); struct2cell(demandATP.SYNPD)], names, 1);
    demandATP=rmfield(demandATP,{'SYN','SYNPD'});
end
types=fieldnames(demandATP);
for i=1:length(types)
    models=fieldnames(demandATP.(types{i}));
    models=models(~contains(models,'constrain'));% ignore unconstrained models
    % 1. plot for two models comparison
    a = figure('units','normalized','outerposition',[0 0 1 1]);
    
    %The reactions and X-axis (ATP demand) are the same for all models
    rxns = {'ATPS4mi', 'glycolysis', 'PGK', 'PYK', 'r0408', ...
        'RE0124C', 'SUCOASm', 'NDPK6', 'NDPK2', 'NMNATr'};
    if any(contains(models,'PD'))
        % PD and control
        if any(contains(models,'1'))
            type.index1=contains(models,'1');
        end
        if any(contains(models,'2'))
            type.index2=contains(models,'2');
        end

        if ~exist('type','var')
            type.index2=true(1,length(models));
        end
        
        for k=1:length(fieldnames(type))
            name=fieldnames(type);
            sets=models(type.(name{k}));
            model1=demandATP.(types{i}).(sets{1});
            model2=demandATP.(types{i}).(sets{2});
            
            customColors1 = [ 0.9255, 0.4784, 0.0314;   % Orange
                0, 0.3725, 0.3765;   % Cyan
                0.9412, 0.6706, 0;   % Yellow
                0, 0.2941, 0.5843;   % deep Blue
                0.2353, 0.2392, 0.6;   % Purple
                0.2196, 0.5059, 0.1843;   % green
                0.416, 0.431, 0.451];   % grey
            
            customColors2 = [0.9569, 0.7137, 0.4706;   % light Orange
                0.6353, 0.8510, 0.8510;   % light cyan
                0.9765, 0.8784, 0.6353;   % light Yellow
                0.5451, 0.7569, 0.9686;   % light blue
                0.6980, 0.6902, 0.9176;   % light Purple
                0.7412, 0.8863, 0.7255; % light green
                0.824, 0.824, 0.824];   % light grey
            
            Pathways={'oxidative phosphorylation (ATPS4mi)', 'glycolysis (PGK & PYK)', ...
                'pentose phosphate pathway (r0408 & RE0124C)', ...
                'citric acid cycle (SUCOASm)',...
                'nucleotide interconversion (NDPK1-10 & UMPK & URIDK3 & ADK1 & r0345 & CYTK1-2)', 'NAD metabolism (NMNATr)',...
                'other bioenergetic reactions'};
            
            % y data for each model
            %             if ~isempty(cell2mat(model1.UMPK))
            ydata = [cell2mat(model1.ATPS4mi(:,1)), cell2mat(model1.glycolysis(:,1)), (cell2mat(model1.r0408(:,1)) ...
                + cell2mat(model1.RE0124C(:,1))), cell2mat(model1.SUCOASm(:,1)),...
                (cell2mat(model1.UMPK(:,1)) + cell2mat(model1.URIDK3(:,1)) + cell2mat(model1.NDPK1(:,1)) + ...
                cell2mat(model1.NDPK2(:,1)) + cell2mat(model1.NDPK3(:,1)) + cell2mat(model1.NDPK4(:,1)) + ...
                cell2mat(model1.NDPK5(:,1)) + cell2mat(model1.NDPK6(:,1)) + cell2mat(model1.NDPK7(:,1)) + ...
                cell2mat(model1.NDPK8(:,1)) + cell2mat(model1.NDPK9(:,1)) + cell2mat(model1.NDPK10(:,1)) + ...
                cell2mat(model1.r0345(:,1)) + cell2mat(model1.ADK1(:,1)) + cell2mat(model1.CYTK1(:,1))+ cell2mat(model1.CYTK2(:,1))),...
                cell2mat(model1.NMNATr(:,1))];
            
            for j = 1:size(ydata,1);
                odata(j,1) = 100 - sum(ydata(j,:));
            end
            ydata(:,7) = odata;
            
            ydataPD= [cell2mat(model2.ATPS4mi(:,1)), cell2mat(model2.glycolysis(:,1)), (cell2mat(model2.r0408(:,1)) ...
                + cell2mat(model2.RE0124C(:,1))), cell2mat(model2.SUCOASm(:,1)),...
                (cell2mat(model2.UMPK(:,1)) + cell2mat(model2.URIDK3(:,1)) + cell2mat(model2.NDPK1(:,1))...
                + cell2mat(model2.NDPK2(:,1))+ cell2mat(model2.NDPK3(:,1)) + cell2mat(model2.NDPK4(:,1))+...
                cell2mat(model2.NDPK5(:,1)) + cell2mat(model2.NDPK6(:,1))+ cell2mat(model2.NDPK7(:,1)) +...
                cell2mat(model2.NDPK8(:,1)) + cell2mat(model2.NDPK9(:,1)) + cell2mat(model2.NDPK10(:,1)) +...
                cell2mat(model2.r0345(:,1)) + cell2mat(model2.ADK1(:,1)) + cell2mat(model2.CYTK1(:,1)) + cell2mat(model2.CYTK2(:,1))),...
                cell2mat(model2.NMNATr(:,1))];
            for m = 1:size(ydataPD,1);
                odataPD(m,1) = 100 - sum(ydataPD(m,:));
            end
            ydataPD(:,7) = odataPD;
            
            %%%%%%
            %%%%%%
            %%%%%% 1 split the stack bars to independent bars
            switch f
                case 'barchart'
                    w1 = 0.4;
                    str = {'%'};
                    col=size(ydata,2);
                    f2= figure('units','normalized','outerposition',[0 0 1 1]);
                    for j = 1:col
                        subplot(col,1,j);
                        % add model1 plot
                        ydata_new=ydata(:,j);
                        bars=bar(d, ydata_new,w1, 'stacked','FaceColor', 'flat');
                        colormap(customColors1);
                        set(bars,'Facecolor',customColors1(j,:));
                        set(gca, 'box','off')
                        for m=1:size(ydata_new,1);
                            xpos = d(m);
                            label=num2str(ydata_new(m),'%.1f');
                            label = strcat(label, '%');
                            hText = text(xpos, ydata_new(m), label);
                            set(hText, 'VerticalAlignment','bottom', 'HorizontalAlignment',...
                                'center','FontSize',10, 'Color','k', 'FontWeight', 'Bold');
                        end
                        
                        hold on
                        % add model2 plot
                        ydataPD_new=ydataPD(:,j);
                        bars=bar(d+4, ydataPD_new,w1, 'stacked','FaceColor', 'flat');
                        colormap(customColors2);
                        set(bars,'Facecolor',customColors2(j,:));
                        ylim([0,max(ydataPD_new)+10]);
                        set(gca, 'box','off')
                        for m=1:size(ydataPD_new,1);
                            xpos = d(m)+4;
                            label=num2str(ydataPD_new(m),'%.1f');
                            label = strcat(label, '%');
                            hText = text(xpos, ydataPD_new(m), label);
                            set(hText, 'VerticalAlignment','bottom', 'HorizontalAlignment',...
                                'center','FontSize',10, 'Color','k', 'FontWeight', 'Bold');
                        end
                        ylim([0,max([ydata_new;ydataPD_new])+5]);
                        xlim([5,110]);
                        set(gca,'XTickLabelRotation',0);
                        title(Pathways{j},'FontSize',12);
                        set(gca,'FontSize',12, 'XTick',...
                            [10 20 30 40 50 60 70 80 90 100], ...
                            'XTickLabel', {...
                            '        C      PD','        C      PD','        C      PD','        C      PD',...
                            '        C      PD','        C      PD','        C      PD','        C      PD',...
                            '        C      PD','        C      PD','        C      PD','        C      PD'});
                        
                        if j~=7
                            set(gca,'XTick',[]);
                        end
                        if j==7
                            xlabel({'Energy demand (ATPM,umol/gDW/h)'},'Position', [55 -5.5 0],'FontSize',20)
                            figure1 = f2;
                            %                         annotate
                            annotate_v1
                        end
                    end
                    % Give common xlabel, ylabel and title to your figure
                    han=axes(f2,'visible','off');
                    han.Title.Visible='on';
                    han.XLabel.Visible='on';
                    han.YLabel.Visible='on';
                    ylabel(han,'ATP Contribution percentage (%)','FontSize',20);
                    sets{2} = strrep(sets{2}, '_', '\_');
                    title(han,[sets{1} ' vs ' sets{2}],'FontSize',20);
                    %%%%%%%%
            end
        end
        
    end
end
%%%%% 2. compare SYN and ASYN
models=fieldnames(demandATP.(types{1}));
models=models(~contains(models,'constrain'));% ignore unconstrained models
for m = 1:length(models)
    sets = fieldnames(demandATP.(types{1}));
    model1name=sets{m};
    model1=demandATP.(types{1}).(sets{m});
    sets = fieldnames(demandATP.(types{2}));
    model2name=sets{m};
    model2=demandATP.(types{2}).(sets{m});
    % y data for each model
    %             if ~isempty(cell2mat(model1.UMPK))
    ydata = [cell2mat(model1.ATPS4mi(:,1)), cell2mat(model1.glycolysis(:,1)), (cell2mat(model1.r0408(:,1)) ...
        + cell2mat(model1.RE0124C(:,1))), cell2mat(model1.SUCOASm(:,1)),...
        (cell2mat(model1.UMPK(:,1)) + cell2mat(model1.URIDK3(:,1)) + cell2mat(model1.NDPK1(:,1)) + ...
        cell2mat(model1.NDPK2(:,1))+ cell2mat(model1.NDPK3(:,1)) + cell2mat(model1.NDPK4(:,1))+ ...
        cell2mat(model1.NDPK5(:,1)) + cell2mat(model1.NDPK6(:,1)) + cell2mat(model1.NDPK7(:,1)) + ...
        cell2mat(model1.NDPK8(:,1)) + cell2mat(model1.NDPK9(:,1)) + cell2mat(model1.NDPK10(:,1)) + ...
        cell2mat(model1.r0345(:,1)) + cell2mat(model1.ADK1(:,1)) + cell2mat(model1.CYTK1(:,1))+ cell2mat(model1.CYTK2(:,1))),...
        cell2mat(model1.NMNATr(:,1))];
    
    for j = 1:size(ydata,1);
        odata(j,1) = 100 - sum(ydata(j,:));
    end
    ydata(:,7) = odata;
    
    ydataPD= [cell2mat(model2.ATPS4mi(:,1)), cell2mat(model2.glycolysis(:,1)), (cell2mat(model2.r0408(:,1)) ...
        + cell2mat(model2.RE0124C(:,1))), cell2mat(model2.SUCOASm(:,1)),...
        (cell2mat(model2.UMPK(:,1)) + cell2mat(model2.URIDK3(:,1)) + cell2mat(model2.NDPK1(:,1))...
        + cell2mat(model2.NDPK2(:,1))+ cell2mat(model2.NDPK3(:,1)) + cell2mat(model2.NDPK4(:,1))+...
        cell2mat(model2.NDPK5(:,1)) + cell2mat(model2.NDPK6(:,1))+ cell2mat(model2.NDPK7(:,1)) +...
        cell2mat(model2.NDPK8(:,1)) + cell2mat(model2.NDPK9(:,1)) + cell2mat(model2.NDPK10(:,1)) +...
        cell2mat(model2.r0345(:,1)) + cell2mat(model2.ADK1(:,1)) + cell2mat(model2.CYTK1(:,1)) + cell2mat(model2.CYTK2(:,1))),...
        cell2mat(model2.NMNATr(:,1))];
    for m = 1:size(ydataPD,1);
        odataPD(m,1) = 100 - sum(ydataPD(m,:));
    end
    ydataPD(:,7) = odataPD;
    
    %%%%%%
    %%%%%%
    %%%%%% 1 split the stack bars to independent bars
    switch f
        case 'barchart'
            w1 = 0.4;
            str = {'%'};
            col=size(ydata,2);
            f2= figure('units','normalized','outerposition',[0 0 1 1]);
            for j = 1:col
                subplot(col,1,j);
                % add model1 plot
                ydata_new=ydata(:,j);
                bars=bar(d, ydata_new,w1, 'stacked','FaceColor', 'flat');
                colormap(customColors1);
                set(bars,'Facecolor',customColors1(j,:));
                set(gca, 'box','off')
                for m=1:size(ydata_new,1);
                    xpos = d(m);
                    label=num2str(ydata_new(m),'%.1f');
                    label = strcat(label, '%');
                    hText = text(xpos, ydata_new(m), label);
                    set(hText, 'VerticalAlignment','bottom', 'HorizontalAlignment',...
                        'center','FontSize',10, 'Color','k', 'FontWeight', 'Bold');
                end
                
                hold on
                % add model2 plot
                ydataPD_new=ydataPD(:,j);
                bars=bar(d+4, ydataPD_new,w1, 'stacked','FaceColor', 'flat');
                colormap(customColors2);
                set(bars,'Facecolor',customColors2(j,:));
                ylim([0,max(ydataPD_new)+10]);
                set(gca, 'box','off')
                for m=1:size(ydataPD_new,1);
                    xpos = d(m)+4;
                    label=num2str(ydataPD_new(m),'%.1f');
                    label = strcat(label, '%');
                    hText = text(xpos, ydataPD_new(m), label);
                    set(hText, 'VerticalAlignment','bottom', 'HorizontalAlignment',...
                        'center','FontSize',10, 'Color','k', 'FontWeight', 'Bold');
                end
                ylim([0,max([ydata_new;ydataPD_new])+5]);
                xlim([5,110]);
                set(gca,'XTickLabelRotation',0);
                title(Pathways{j},'FontSize',12);
                set(gca,'FontSize',12, 'XTick',...
                    [10 20 30 40 50 60 70 80 90 100], ...
                    'XTickLabel', {...
                    '        SYN     ASYN','        SYN     ASYN','        SYN     ASYN','        SYN     ASYN',...
                    '        SYN     ASYN','        SYN     ASYN','        SYN     ASYN','        SYN     ASYN',...
                    '        SYN     ASYN','        SYN     ASYN','        SYN     ASYN','        SYN     ASYN'});
                
                if j~=7
                    set(gca,'XTick',[]);
                end
                if j==7
                    xlabel({'Energy demand (ATPM,umol/gDW/h)'},'Position', [55 -5.5 0],'FontSize',20)
                    figure1 = f2;
                    %                         annotate
                    annotate_v1
                end
            end
            % Give common xlabel, ylabel and title to your figure
            han=axes(f2,'visible','off');
            han.Title.Visible='on';
            han.XLabel.Visible='on';
            han.YLabel.Visible='on';
            ylabel(han,'ATP Contribution percentage (%)','FontSize',20);
            title(han,[model1name ' vs ' model2name],'FontSize',20);
            %%%%%%%%
    end
    
end

end





