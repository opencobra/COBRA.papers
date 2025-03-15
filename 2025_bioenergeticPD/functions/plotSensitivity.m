function plotSensitivity(Allmodels,ignoreUnconstrain,method, ATPMvalue)
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
    if any(contains(modelnum,'PD'))
        % PD and control
        if length(modelnum) == 2
            % control and PD
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %      oxygen sensitivity
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            xlabel('The uptake of Oxygen (\mumol/gDW/h)', 'FontSize', 14);
            ylabel('The flux through ATP synthase (\mumol/gDW/h)','FontSize', 14);
            title([modelnum{1} ' vs ' modelnum{2}],'FontSize', 14)
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
            legend(modelnum{1},['change point of ', modelnum{1}],modelnum{2},['change point of ', modelnum{2}],'box','off','FontSize', 14, 'Location', 'best');
            clear m n OCR ATPS4mi_ox
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %    glucose sensitivity
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % control and PD
            model1=Allmodels.(types{i}).(modelnum{1});
            model2=Allmodels.(types{i}).(modelnum{2});
            % change the ATPM value
            rxn=model1.rxns(contains(model1.rxns,'ATPM'));
            model1 = changeRxnBounds(model1,rxn,ATPMvalue,'l');
            %
            GCR = linspace(model1.lb(find(ismember(model1.rxns, 'EX_glc_D[e]'))),0,10);
            for m =1:length(GCR);
                model1 = changeRxnBounds(model1,'EX_glc_D[e]',GCR(m),'b');
                model1 = changeObjective(model1,'ATPM');
                switch method
                    case 'FBA'
                        FBAenergetics = optimizeCbModel(model1,'max', 1e-6);
                        ATPS4mi_ox(m,1) = FBAenergetics.x(find(contains(model1.rxns, 'ATPS4mi')));
                    case 'eFBA'
                        try
                            FBAenergetics=entropicFluxBalanceAnalysis(model1,param);
                            ATPS4mi_ox(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'ATPS4mi')));
                            lac(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'EX_lac_L[e]')));
                            gln(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'EX_gln_L[e]')));
                            gal(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'EX_gal[e]')));
                        catch
                            ATPS4mi_ox(m,1)=0;
                            lac(m,1)=0;
                            gln(m,1)=0;
                            gal(m,1)=0;
                        end
                end
            end
            % Define a fashionable color palette
            customColors1 = [ 0.9255, 0.4784, 0.0314;   % Orange
                0, 0.2941, 0.5843;   % deep Blue
                0, 0.3725, 0.3765;   % Cyan
                0.9412, 0.6706, 0];   % Yellow
            
            
            customColors2 = [0.9569, 0.7137, 0.4706;   % light Orange
                0.5451, 0.7569, 0.9686;   % light blue
                0.6353, 0.8510, 0.8510;   % light cyan
                0.9765, 0.8784, 0.6353];   % light Yellow
            
            
            subplot(2,2,i+2);
            plot(GCR, ATPS4mi_ox, 'LineWidth', 2, 'Color',customColors1(1,:))
            % add another model
            % change the ATPM value
            rxn=model2.rxns(contains(model2.rxns,'ATPM'));
            model2 = changeRxnBounds(model2,rxn,ATPMvalue,'l');
            %
            GCR = linspace(model2.lb(find(ismember(model2.rxns, 'EX_glc_D[e]'))),0,10);
            for n =1:length(GCR);
                model2 = changeRxnBounds(model2,'EX_glc_D[e]',GCR(n),'b');
                model2 = changeObjective(model2,'ATPM');
                
                switch method
                    case 'FBA'
                        FBAenergetics = optimizeCbModel(model2,'max', 1e-6);
                        ATPS4mi_ox(n,1) = FBAenergetics.x(find(contains(model2.rxns, 'ATPS4mi')));
                    case 'eFBA'
                        try
                            FBAenergetics=entropicFluxBalanceAnalysis(model2,param);
                            ATPS4mi_ox(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'ATPS4mi')));
                            lac(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'EX_lac_L[e]')));
                            gln(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'EX_gln_L[e]')));
                            gla(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'EX_gal[e]')));
                        catch
                            ATPS4mi_ox(n,1)=0;
                            lac(n,1)=0;
                            gln(n,1)=0;
                            ga(n,1)=0;
                        end
                end
            end
            hold on
            plot(GCR, ATPS4mi_ox, 'LineWidth', 2, 'Color',customColors1(2,:))
            xlabel('The uptake of Glucose (\mumol/gDW/h)', 'FontSize', 14);
            ylabel('The flux through ATP synthase (\mumol/gDW/h)','FontSize', 14);
            ylim([min(ATPS4mi_ox)-0.5,max(ATPS4mi_ox)+1])
            title([modelnum{1} ' vs ' modelnum{2}],'FontSize', 14)
            legend(modelnum{1},modelnum{2},'box','off','FontSize', 16, 'Location', 'best')
            clear m n OCR ATPS4mi_ox
        end
        
    end
end


