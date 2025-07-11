function plotChangingTrends(demandATP,f,d)
%   robustness analysis
%   drawing the flux changing trends with different ATPM value
if ~exist('f','var')
    f='barchart';
end

rxns = {'ATPS4mi', 'glycolysis', 'PGK', 'PYK', 'r0408', ...
    'RE0124C', 'SUCOASm', 'NDPK6', 'NDPK2', 'NMNATr'};
customColors1 = [ 0.9255, 0.4784, 0.0314;   % Orange
    0, 0.3725, 0.3765;   % Cyan
    0.9412, 0.6706, 0;   % Yellow
    0, 0.2941, 0.5843;   % deep Blue
    0.2353, 0.2392, 0.6;   % Purple
    0.2196, 0.5059, 0.1843;   % green
    0.416, 0.431, 0.451];   % grey

Pathways={'oxidative phosphorylation (ATPS4mi)', 'glycolysis (PGK & PYK)', ...
                'pentose phosphate pathway (r0408 & RE0124C)', ...
                'citric acid cycle (SUCOASm)',...
                'nucleotide interconversion (NDPK1-10 & UMPK & URIDK3 & ADK1 & r0345 & CYTK1-2)', 'NAD metabolism (NMNATr)',...
                'other bioenergetic reactions'};
clear ydata ydata_new
%%%%%%%%%%
types=fieldnames(demandATP);
for i=1:length(types)
    models=fieldnames(demandATP.(types{i}));
    models=models(~contains(models,'constrain'));% ignore unconstrained models
    for j=1:length(models)
        if any(d)
            d = d;
        else
            d=10:10:100;
        end
        %The reactions and X-axis (ATP demand) are the same for all models
        model=demandATP.(types{i}).(models{j});
        
       
        % y data (flux value of each pathway) for each model
        ydata = [cell2mat(model.ATPS4mi(:,2)), cell2mat(model.PGK(:,2))+ cell2mat(model.PYK(:,2)), (cell2mat(model.r0408(:,2)) ...
            + cell2mat(model.RE0124C(:,2))), cell2mat(model.SUCOASm(:,2)),...
            (cell2mat(model.UMPK(:,2)) + cell2mat(model.URIDK3(:,2)) + cell2mat(model.NDPK1(:,2)) + cell2mat(model.NDPK2(:,2))+ cell2mat(model.NDPK3) + cell2mat(model.NDPK4(:,2))+ cell2mat(model.NDPK5(:,2)) + cell2mat(model.NDPK6) + cell2mat(model.NDPK7) + cell2mat(model.NDPK8) + cell2mat(model.NDPK9(:,2)) + cell2mat(model.NDPK10) + cell2mat(model.r0345(:,2)) + cell2mat(model.ADK1(:,2))),...
            cell2mat(model.NMNATr(:,2))];
        
        switch f
            case 'line'
                col=size(Pathways,2);  
                subplot(2,2,i)
                for k = 1:col
                    ydata_new=ydata(:,k);
                    % Plot each line
                    plot(d, ydata_new, '-o', 'DisplayName', ['Line ' num2str(k)],'Color', customColors1(k,:),'LineWidth', 1.5);
                    % Add label at the last data point
                    shortPath= regexp(Pathways{k}, '^[^()]+', 'match', 'once');
                    text(d(end), ydata_new(end), shortPath, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','Color', customColors1(k,:)); %, 'Color', customColors1(k,:)
                    hold on
                end
                % Add labels and title
                xlabel('Energy demand (ATPM, umol/h/gDW)','FontSize',14);
                ylabel('Flux value (umol/h/gDW)','FontSize',14);
                title(['Flux changes on different pathways in ', models{j}],'FontSize', 14);
                
                %
                changePoint = 60;
                plot([changePoint,changePoint],[0, 90], 'r--', 'LineWidth', 1.5);
%                 text(changePoint,100, ['X=' num2str(changePoint)], ...
%                     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                if i==1
                legend([Pathways 'changeLine'], 'Location', 'best','Box','off','FontSize', 10);
                end
        end
    end
end

