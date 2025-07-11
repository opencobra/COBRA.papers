function plotRespirationRobustness(Allmodels,ignoreUnconstrain,method, ATPMvalue)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~isstruct(Allmodels)
    error('check input models')
end

if ~exist('ignoreUnconstrain','var')
    ignoreUnconstrain=0;
end

if ~exist('method','var')
    method='FBA';
else
    param.debug = 1;
    param.feasTol=1e-7;
end

if ~exist('ATPMvalue','var')
    error('please set ATP requirement value')
end

if any(contains(fieldnames(Allmodels),{'SYNPD'}))
    names = [fieldnames(Allmodels.SYN); fieldnames(Allmodels.SYNPD)];
    Allmodels.allSYN = cell2struct([struct2cell(Allmodels.SYN); struct2cell(Allmodels.SYNPD)], names, 1);
    names = [fieldnames(Allmodels.ASYN); fieldnames(Allmodels.ASYNPD)];
    Allmodels.allASYN = cell2struct([struct2cell(Allmodels.ASYN); struct2cell(Allmodels.ASYNPD)], names, 1);
    
    Allmodels=rmfield(Allmodels,{'SYN','SYNPD','ASYN','ASYNPD'});
end

types=fieldnames(Allmodels);
for i=1:length(types)
    models=Allmodels.(types{i});
    modelnum=fieldnames(models);
    if ignoreUnconstrain ~= 0
        modelnum=modelnum(~contains(modelnum,'constrain'));% ignore unconstrained models
    end
    if ~any(contains(modelnum,'ASTRO'))
        if any(contains(modelnum,'PD'))
            % PD and control
            if length(modelnum) == 4
                group.index1=contains(modelnum,'1');
                group.index2=contains(modelnum,'2');
                d = figure('units','normalized','outerposition',[0 0 1 1]);
                for k=1:length(fieldnames(group))
                    name=fieldnames(group);
                    sets=modelnum(group.(name{k}));
                    model1=Allmodels.(types{i}).(sets{1});
                    model2=Allmodels.(types{i}).(sets{2});
                    OCR = [model1.lb(find(ismember(model1.rxns, 'EX_o2[e]'))):...
                        0]';%model1.ub(find(ismember(model1.rxns, 'EX_o2[e]')))
                    for m =1:length(OCR);
                        model1 = changeRxnBounds(model1,'EX_o2[e]',OCR(m,1),'b');
                        model1 = changeObjective(model1,'ATPM');
                        switch method
                            case 'FBA'
                                FBAenergetics = optimizeCbModel(model1,'max', 1e-6);
                                %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(m,1) = FBAenergetics.x(find(contains(model1.rxns, 'ATPS4mi')));
                                %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                            case 'eFBA'
                                try
                                    FBAenergetics=entropicFluxBalanceAnalysis(model1,param);
                                    %                             CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                    ATPS4mi_ox(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'ATPS4mi')));
                                catch
                                    %                             CYOOm3i_ox(k,1)=0;
                                    ATPS4mi_ox(m,1)=0;
                                end
                        end
                    end
                    subplot(1,2,k);
                    plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                    hold on
                    clear m n ATPS4mi_ox
                    %%%%%%%%%%%%
                    % add another model
                    for n =1:length(OCR);
                        model2 = changeRxnBounds(model2,'EX_o2[e]',OCR(n,1),'b');
                        model2 = changeObjective(model2,'ATPM');
                        switch method
                            case 'FBA'
                                FBAenergetics = optimizeCbModel(model2,'max', 1e-6);
                                %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(n,1) = FBAenergetics.x(find(contains(model2.rxns, 'ATPS4mi')));
                                %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                            case 'eFBA'
                                try
                                    FBAenergetics=entropicFluxBalanceAnalysis(model2,param);
                                    %                             CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                    ATPS4mi_ox(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'ATPS4mi')));
                                catch
                                    %                             CYOOm3i_ox(k,1)=0;
                                    ATPS4mi_ox(n,1)=0;
                                end
                        end
                    end
                    plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                    xlabel('Oxygen exchange flux (\mumol/gDW/h)', 'FontSize', 24);
                    ylabel('ATP synthase flux (\mumol/gDW/h)','FontSize', 24);
                    title([sets{1} ' vs ' sets{2}])
                    legend(sets{1},sets{2},'box','off','FontSize', 14)
                    clear m n OCR ATPS4mi_ox
                end
                clear group d name sets model1 model2
                %
                %
                % old and new
                group.index1=~contains(modelnum,'PD');
                group.index2=contains(modelnum,'PD');
                d = figure('units','normalized','outerposition',[0 0 1 1]);
                for k=1:length(fieldnames(group))
                    name=fieldnames(group);
                    sets=modelnum(group.(name{k}));
                    model1=Allmodels.(types{i}).(sets{1});
                    model2=Allmodels.(types{i}).(sets{2});
                    OCR = [model1.lb(find(ismember(model1.rxns, 'EX_o2[e]'))):...
                        0]';%model1.ub(find(ismember(model1.rxns, 'EX_o2[e]')))
                    for m =1:length(OCR);
                        model1 = changeRxnBounds(model1,'EX_o2[e]',OCR(m,1),'b');
                        model1 = changeObjective(model1,'ATPM');
                        switch method
                            case 'FBA'
                                FBAenergetics = optimizeCbModel(model1,'max', 1e-6);
                                %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(m,1) = FBAenergetics.x(find(contains(model1.rxns, 'ATPS4mi')));
                                %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                            case 'eFBA'
                                try
                                    FBAenergetics=entropicFluxBalanceAnalysis(model1,param);
                                    %                             CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                    ATPS4mi_ox(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'ATPS4mi')));
                                catch
                                    %                             CYOOm3i_ox(k,1)=0;
                                    ATPS4mi_ox(m,1)=0;
                                end
                        end
                    end
                    
                    subplot(1,2,k);
                    plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                    hold on
                    %%%%%%%%%%%%
                    % add another model
                    OCR = [model2.lb(find(ismember(model2.rxns, 'EX_o2[e]'))):...
                        0]';%model2.ub(find(ismember(model2.rxns, 'EX_o2[e]')))
                    for n =1:length(OCR);
                        model2 = changeRxnBounds(model2,'EX_o2[e]',OCR(n,1),'b');
                        model2 = changeObjective(model2,'ATPM');
                        switch method
                            case 'FBA'
                                FBAenergetics = optimizeCbModel(model2,'max', 1e-6);
                                %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(n,1) = FBAenergetics.x(find(contains(model2.rxns, 'ATPS4mi')));
                                %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                            case 'eFBA'
                                try
                                    FBAenergetics=entropicFluxBalanceAnalysis(model2,param);
                                    %                             CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                    ATPS4mi_ox(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'ATPS4mi')));
                                catch
                                    %                             CYOOm3i_ox(k,1)=0;
                                    ATPS4mi_ox(n,1)=0;
                                end
                                
                        end
                    end
                    plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                    xlabel('Oxygen exchange flux (\mumol/gDW/h)', 'FontSize', 24);
                    ylabel('ATP synthase flux (\mumol/gDW/h)','FontSize', 24);
                    title([sets{1} ' vs ' sets{2}])
                    legend(sets{1},sets{2},'box','off','FontSize', 14)
                    clear m n OCR ATPS4mi_ox
                end       
            elseif length(modelnum) == 2
                % control and PD
                model1=Allmodels.(types{i}).(modelnum{1});
                model2=Allmodels.(types{i}).(modelnum{2});
                % change the ATPM value
                rxn=model1.rxns(contains(model1.rxns,'ATPM'));
                model1 = changeRxnBounds(model1,rxn,ATPMvalue,'l');
                OCR = [model1.lb(find(ismember(model1.rxns, 'EX_o2[e]'))):...
                    0]';%model1.ub(find(ismember(model1.rxns, 'EX_o2[e]')))
                for m =1:length(OCR);
                    model1 = changeRxnBounds(model1,'EX_o2[e]',OCR(m,1),'b');
                    model1 = changeObjective(model1,'ATPM');
                    switch method
                        case 'FBA'
                            FBAenergetics = optimizeCbModel(model1,'max', 1e-6);
                            ATPS4mi_ox(m,1) = FBAenergetics.x(find(contains(model1.rxns, 'ATPS4mi')));
                        case 'eFBA'
                            try
                                FBAenergetics=entropicFluxBalanceAnalysis(model1,param);
                                ATPS4mi_ox(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'ATPS4mi')));
                            catch
                                ATPS4mi_ox(m,1)=0;
                            end
                    end
                end
                
                % Define a fashionable color palette
                colors = [0, 0.2941, 0.5843;  %  blue
                    0.9255, 0.4784, 0.0314]; % orange
                
                subplot(2,2,i);
                % Create a bar plot with a custom color scheme
                p1=plot(OCR, ATPS4mi_ox, 'LineWidth', 2);
                p1.Color = colors(1,:);
                % add changing point
                diff_data = diff(ATPS4mi_ox);
                [~, idx] = max(abs(diff_data));
                changePoint = OCR(idx-2);
                specificY = ATPS4mi_ox(idx-2);
                hold on
                if ~isempty(changePoint)
                    plot([changePoint,changePoint],[0, specificY], 'r--', 'LineWidth', 1.5);
                    text(changePoint,specificY, ['O2uptake =' num2str(changePoint)], ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                end
                hold on
                clear OCR
                %%%%%%%%%%%%
                % add another model
                % change the ATPM value
                rxn=model2.rxns(contains(model2.rxns,'ATPM'));
                model2 = changeRxnBounds(model2,rxn,ATPMvalue,'l');
                %
                OCR = [model2.lb(find(ismember(model2.rxns, 'EX_o2[e]'))):...
                    0]';%model1.ub(find(ismember(model1.rxns, 'EX_o2[e]')))
                for n =1:length(OCR);
                    model2 = changeRxnBounds(model2,'EX_o2[e]',OCR(n,1),'b');
                    model2 = changeObjective(model2,'ATPM');
                    switch method
                        case 'FBA'
                            FBAenergetics = optimizeCbModel(model2,'max', 1e-6);
                            ATPS4mi_ox(n,1) = FBAenergetics.x(find(contains(model2.rxns, 'ATPS4mi')));
                        case 'eFBA'
                            try
                                FBAenergetics=entropicFluxBalanceAnalysis(model2,param);
                                ATPS4mi_ox(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'ATPS4mi')));
                            catch
                                ATPS4mi_ox(n,1)=0;
                            end
                            
                    end
                end
                % Create a bar plot with a custom color scheme
                p1=plot(OCR, ATPS4mi_ox, 'LineWidth', 2);
                p1.Color = colors(2,:);
                xlabel('Oxygen exchange flux (\mumol/gDW/h)', 'FontSize', 14);
                ylabel('ATP synthase flux (\mumol/gDW/h)','FontSize', 14);
                title([modelnum{1} ' vs ' modelnum{2}],'FontSize', 16)
                % add changing point
                % add changing point
                diff_data = diff(ATPS4mi_ox);
                [~, idx] = max(abs(diff_data));
                changePoint = OCR(idx-2);
                specificY = ATPS4mi_ox(idx-2);
                hold on
                if ~isempty(changePoint)
                    plot([changePoint,changePoint],[0, specificY], 'r--', 'LineWidth', 1.5);
                    text(changePoint,specificY, ['O2 uptake =' num2str(changePoint)], ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                end
                %
                legend(modelnum{1},['changepoint of ', modelnum{1}],modelnum{2},['changepoint of ', modelnum{2}],'box','off','FontSize', 14, 'Location', 'best');
                clear m n OCR ATPS4mi_ox
            end
        else
            if ~any(contains(modelnum,'constrain'))
                %old model
                model1 = Allmodels.(types{i}).(modelnum{1});
                % Define the flux ranges based on the lower and upper bounds of the
                % corresponding models
                % Using same bounds for control and PD models.
                OCR = [model1.lb(find(ismember(model1.rxns, 'EX_o2[e]'))):...
                    model1.ub(find(ismember(model1.rxns, 'EX_o2[e]')))]';
                for k =1:length(OCR);
                    model1 = changeRxnBounds(model1,'EX_o2[e]',OCR(i,1),'b');
                    model1 = changeObjective(model1,'ATPM');
                    switch method
                        case 'FBA'
                            FBAenergetics = optimizeCbModel(model1,'max', 1e-6);
                            %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                            ATPS4mi_ox(k,1) = FBAenergetics.x(find(contains(model1.rxns, 'ATPS4mi')));
                            %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                        case 'eFBA'
                            try
                                FBAenergetics=entropicFluxBalanceAnalysis(model1,param);
                                ATPS4mi_ox(k,1) = FBAenergetics.v(find(contains(model1.rxns, 'ATPS4mi')));
                            catch
                                ATPS4mi_ox(k,1)=0;
                            end
                    end
                end
                d = figure('units','normalized','outerposition',[0 0 1 1]);
                plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                hold on
                %
                %
                %new model
                model2 = Allmodels.(types{i}).(modelnum{2});
                % Define the flux ranges based on the lower and upper bounds of the
                % corresponding models
                % Using same bounds for control and PD models.
                clear ATPS4mi_ox
                for k =1:length(OCR);
                    model2 = changeRxnBounds(model2,'EX_o2[e]',OCR(k,1),'b');
                    model2 = changeObjective(model2,'ATPM');
                    switch method
                        case 'FBA'
                            FBAenergetics = optimizeCbModel(model2,'max', 1e-6);
                            if FBAenergetics.stat ==1
                                %CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(k,1) = FBAenergetics.x(find(contains(model2.rxns, 'ATPS4mi')));
                                %ATPM_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                            else
                                ATPS4mi_ox(k,1) =0;
                            end
                        case 'eFBA'
                            try
                                FBAenergetics=entropicFluxBalanceAnalysis(model2,param);
                                %CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(k,1) = FBAenergetics.v(find(contains(model2.rxns, 'ATPS4mi')));
                            catch
                                %CYOOm3i_ox(k,1)=0;
                                ATPS4mi_ox(k,1)=0;
                            end
                    end
                end
                plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                xlabel('Oxygen exchange flux (\mumol/gDW/h)', 'FontSize', 24);
                ylabel('ATP synthase flux (\mumol/gDW/h)','FontSize', 24);
                hold off
                title([modelnum{1} ' vs ' modelnum{2}])
                legend(modelnum{1},modelnum{2},'box','off','FontSize', 14)
                ylim(gca,[2, 7])
            end
        end
    end
end

