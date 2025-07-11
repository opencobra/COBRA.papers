function plotExFluxes(Allmodels,obj,method)
%plot the flux value for exchange rxns
modeltypes=fieldnames(Allmodels);

for i=1:length(modeltypes)
    models=Allmodels.(modeltypes{i});
    modelnum=fieldnames(models);
    modelnum=modelnum(~contains(modelnum,'constrain'));% ignore constrained models
    for j=1:length(modelnum)
        model=Allmodels.(modeltypes{i}).(modelnum{j});
        if ~isfield(model,'ExchRxnBool')
            model=findSExRxnInd(model);
        end
        rxn=model.rxns(contains(model.rxns,obj));
        model = changeObjective(model,rxn);
        switch method
            case 'FBA'
                solBMa =optimizeCbModel(model,'max',1e-4);
                solBMa.obj = solBMa.f; % purpose of renaming fields?
                solBMa.full = solBMa.x;
                %setting the fluxes below eucNorm to zero
                for k=1:length(solBMa.full)
                    if abs(solBMa.full(k))< eucNorm % threshold applied to solBMa.x but flux splits computed with solBMa.full
                        solBMa.full(k)=0;
                    end
                end
            case {'eFBA','entropicFBA'}
                clear param
                param.solver = 'mosek';
                param.internalNetFluxBounds='original';
                param.externalNetFluxBounds='original';
                param.method='fluxes';
                param.printLevel=2;
                param.debug = 1;
                param.feasTol=1e-7;
                [solBMa,~] = entropicFluxBalanceAnalysis(model,param);
                if solBMa.stat==1
                    solBMa.full = solBMa.v;
                    for k=1:length(solBMa.full)
                        if abs(solBMa.full(k))< 1e-4
                            solBMa.full(k)=0;
                        end
                    end
                else
                    error(' EPproblem is not feasible')
                end
        end
%         AllEX=model.rxns(model.ExchRxnBool);
        ind1 = (model.ExchRxnBool & (solBMa.full<0));
        RealUptake = [model.rxns(ind1), table(solBMa.full(ind1))];
        RealUptake.Properties.VariableNames = {'RxnID','fluxValue'};
        ind2 = (model.ExchRxnBool & (solBMa.full>0));
        RealExcrete = [model.rxns(ind2), table(solBMa.full(ind2))];
        RealExcrete.Properties.VariableNames = {'RxnID','fluxValue'};
        
        % ranged dot plot for flux value of uptake and Excrete rxns
        % d = figure;
        %% uptake flux
        subplot(2,2,i)
        [A,B]=ismember(model.rxns,RealUptake.RxnID);
        RealUptake.RxnName(B(A))=model.rxnNames(A);
        labels = RealUptake.RxnName; % Labels for the y-axis
        central_values = RealUptake.fluxValue; % Central calculated values 
        % Number of data points
        n = length(labels);
        hold on;
        % Loop through each data point to plot the central values and ranges (lb, ub)
        % uptake flux <0
        minFlux=min(RealUptake.fluxValue);
        for m = 1:n
            % Plot the central value as a dot
            plot(central_values(m), m, 'bo', 'MarkerFaceColor', 'b');
            % Plot the range as a line
            lb = model.lb(ismember(model.rxns,RealUptake.RxnID{m}));
            ub = model.ub(ismember(model.rxns,RealUptake.RxnID{m}));
            if lb <= minFlux
                startpoint = minFlux - 5;
            else
                startpoint = lb;
            end
            if ub >= abs(minFlux)
                endpoint = abs(minFlux) + 5;
            else
                endpoint = ub;
            end
%             line([startpoint, endpoint], [m, m], 'Color', 'r', 'LineWidth', 2);
            line([startpoint, endpoint], [m, m], 'Color', 'black', 'LineWidth', 2);  
            % Add labels for the dots
            dot_label = char(string(central_values(m))); 
            text(central_values(m) + 0.4, m +0.4, dot_label, 'FontSize', 6, 'FontName', 'Arial', 'VerticalAlignment', 'middle');
        end
       
        % Customize the plot       
        set(gca, 'YTick', 1:n, 'YTickLabel', strrep(labels, '_', '\_'), 'FontSize', 8,'FontName', 'Arial'); % Set y-axis labels, specify '_';        
        xlabel('fluxValue (umol/gDW/h)','FontSize',14);
        ylabel('Uptake reactions','FontSize',14);
        title(['Ranged Dot Plot',' of ',strrep(modelnum{j}, '_', '\_')]);
        grid off;
        ylim([0,n+1])
        xlim([minFlux-10, abs(minFlux)+10]); % Set x-axis limits
        % add info line
        plot([0,0],[0, n], 'r--', 'LineWidth', 1.5);
        
        f=round(minFlux - 5,0);
        plot([f,f],[0, n], 'r--', 'LineWidth', 0.5)
        x_text = f-2; % x-coordinate where you want the text
        y_text = -2; % y-coordinate where you want the text (e.g., middle of the line)
        % Add text on the line
        text(x_text, y_text, ['lb < ', char(string(f))], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'FontWeight', 'bold', 'FontName', 'Arial');
        
        f=round(abs(minFlux) + 5,0);
        plot([f,f],[0, n], 'r--', 'LineWidth', 0.5)
        x_text = f-2; % x-coordinate where you want the text
        y_text = -2; % y-coordinate where you want the text (e.g., middle of the line)
        % Add text on the line
        text(x_text, y_text, ['ub > ', char(string(f))], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'FontWeight', 'bold', 'FontName', 'Arial');
        
        % Remove whitespace between subplots
        set(gcf, 'PaperPositionMode', 'auto');

        if 0
        %% excrete flux
        subplot(1,2,2)
        [A,B]=ismember(model.rxns,RealExcrete.RxnID);
        RealExcrete.RxnName(B(A))=model.rxnNames(A);
        labels = RealExcrete.RxnName; % Labels for the y-axis
        central_values = RealExcrete.fluxValue; % Central calculated values, 
        % Number of data points
        n = length(labels);
        hold on;
        % Loop through each data point to plot the central values and ranges (lb, ub)
        % excrete flux > 0
        maxFlux=max(RealExcrete.fluxValue);
        for m = 1:n
            % Plot the central value as a dot
            plot(central_values(m), m, 'bo', 'MarkerFaceColor', 'b');
            % Plot the range as a line
            lb = model.lb(ismember(model.rxns,RealExcrete.RxnID{m}));
            ub = model.ub(ismember(model.rxns,RealExcrete.RxnID{m}));
            if lb <= -maxFlux
                startpoint = -(maxFlux)- 5;
            else
                startpoint = lb;
            end
            if ub >= maxFlux
                endpoint = maxFlux+5;
            else
                endpoint = ub;
            end
%             line([startpoint, endpoint], [m, m], 'Color', 'r', 'LineWidth', 2);
            line([startpoint, endpoint], [m, m], 'Color', 'black', 'LineWidth', 2);
            % Add labels for the dots
            dot_label = char(string(central_values(m))); 
            text(central_values(m) + 0.4, m+0.4, dot_label, 'FontSize', 8, 'FontName', 'Arial', 'VerticalAlignment', 'middle');
        end

        % Customize the plot       
        set(gca, 'YTick', 1:n, 'YTickLabel', strrep(labels, '_', '\_'), 'FontSize', 8,'FontName', 'Arial'); % Set y-axis labels, specify '_';        
        xlabel('fluxValue (umol/gDW/h)','FontSize',14);
        ylabel('Excrete reactions','FontSize',14');
        title(['Ranged Dot Plot',' of ',strrep(modelnum{j}, '_', '\_'), ' on excrete']);
        grid off;
        ylim([0,n+1])
        xlim([-maxFlux-10, maxFlux+10]); % Set x-axis limits
        % add info line
        plot([0,0],[0, n], 'r--', 'LineWidth', 1.5);
        % add note1
        f=round(-maxFlux - 5,0);
        plot([f,f],[0, n], 'r--', 'LineWidth', 0.5)
        x_text = f-2; % x-coordinate where you want the text
        y_text = -4; % y-coordinate where you want the text (e.g., middle of the line)
        % Add text on the line
        text(x_text, y_text, ['lb < ', char(string(f))], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial');
        % add note2
        f=round(maxFlux + 5,0);
        plot([f,f],[0, n], 'r--', 'LineWidth', 0.5)
        x_text = f-2; % x-coordinate where you want the text
        y_text = -4; % y-coordinate where you want the text (e.g., middle of the line)
        % Add text on the line
        text(x_text, y_text, ['ub > ', char(string(f))], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial');
    
        end
    end
end
end

