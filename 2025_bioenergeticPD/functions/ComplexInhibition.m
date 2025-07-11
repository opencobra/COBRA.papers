% Complex I inhibition
All=Allmodels;
% reshape
if any(contains(fieldnames(All),{'SYNPD'}))
    names = [fieldnames(All.SYN); fieldnames(All.SYNPD)];
    All.allSYN = cell2struct([struct2cell(All.SYN); struct2cell(All.SYNPD)], names, 1);
    names = [fieldnames(All.ASYN); fieldnames(All.ASYNPD)];
    All.allASYN = cell2struct([struct2cell(All.ASYN); struct2cell(All.ASYNPD)], names, 1);
    All=rmfield(All,{'SYN','SYNPD','ASYN','ASYNPD'});
end

if ~exist('method','var')
    method='FBA';
else
    param.debug = 1;
    param.feasTol=1e-4;
end

%% explore the complex I inhibition based on default ATPM value (10.6 umol/h/gDW)
types=fieldnames(All);
d = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(types)
    models=All.(types{i});
    modelnum=fieldnames(models);
    % only perform complex I inhibition on unconstrained models
    modelnum=modelnum(~contains(modelnum,'constrain'));% only use unconstrained models!!!
    if ~any(contains(modelnum,'ASTRO'))
        if any(contains(modelnum,'PD'))
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %      PD vs control
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(modelnum) == 4
                group.index1=contains(modelnum,'1');
                group.index2=contains(modelnum,'2');
                if any(contains(types{i},'ASYN'))
                    ATPMvalue = ATPMvalue_ASYN;
                    %                 O2value = O2value_ASYN;
                    %                 ComplexIII= ComplexIII_ASYN;
                else
                    ATPMvalue = ATPMvalue_SYN;
                    %                 O2value = O2value_SYN;
                    %                 ComplexIII= ComplexIII_SYN;
                end
                for k=1:length(fieldnames(group))
                    name=fieldnames(group);
                    sets=modelnum(group.(name{k}));
                    model1=All.(types{i}).(sets{1});
                    model2=All.(types{i}).(sets{2});
                    % control model
                    AtpRates1_ComplexI = zeros(101,1);
                    cIactivityoptimal1 = zeros(101,1);
                    complexIactivity=model1.ub(ismember(model1.rxns,'NADH2_u10mi'));
                    % change ATP requirement
                    rxn3=model1.rxns(contains(model1.rxns,'ATPM'));
                    model1 = changeRxnBounds(model1,rxn3,ATPMvalue,'l');
                    % only set minimal oxygen consumption for control models
                    %                 rxn4=model1.rxns(contains(model1.rxns,'EX_o2[e]'));
                    %                 model1 = changeRxnBounds(model1,rxn4,O2value,'u');
                    % add complex III value???
                    %                 rxn5=model1.rxns(contains(model1.rxns,'CYOR_u10mi'));
                    %                 if any(contains(types{i},'ASYN'))
                    %                     model1 = changeRxnBounds(model1,rxn5,ComplexIII,'u');
                    %                 else
                    %                     model1 = changeRxnBounds(model1,rxn5,ComplexIII,'l');
                    %                 end
                    for m = 0:100;
                        ComplexImodel = model1;
                        % change complex I activity
                        rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
                        ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((m/100)*complexIactivity)),'u');
                        % change objective function
                        rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
                        ComplexImodel = changeObjective(ComplexImodel,rxn2);
                        
                        switch method
                            case 'FBA'
                                FBA_ComplexI = optimizeCbModel(ComplexImodel,'max');
                                FBA_ComplexI.l=FBA_ComplexI.x(find(ismember(ComplexImodel.rxns, rxn1)));
                            case 'eFBA'
                                try
                                    FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);
                                    FBA_ComplexI.f=FBA_ComplexI.v(contains(ComplexImodel.rxns,'ATPS4mi'));
                                    FBA_ComplexI.l=FBA_ComplexI.v(find(ismember(ComplexImodel.rxns, rxn1)));
                                catch
                                    FBA_ComplexI.f=0;
                                    FBA_ComplexI.l=0;
                                end
                        end
                        AtpRates1_ComplexI(m+1) = FBA_ComplexI.f;
                        cIactivityoptimal1(m+1) = FBA_ComplexI.l;
                    end
                    %Plot flux through ATP synthase
                    subplot(1,2,k)
                    ratio=(AtpRates1_ComplexI/(max(AtpRates1_ComplexI)))*100;
                    ratio = round(ratio,2);
                    plot (1:101,ratio, 'LineWidth',2, 'Color',[1 0 1]);
                    
                    % add anotation
                    hold on;  % Retain the current plot
                    changePoint1=min(find(ratio<99.99)-1);
                    specificY1 = ratio(changePoint1);
                    %                 plot(changePoint,specificY,'ro', 'MarkerSize', 10);
                    if ~isempty(changePoint1)
                        plot([changePoint1,changePoint1],[0, specificY1], 'r--', 'LineWidth', 1.5);
                        text(changePoint1,specificY1, ['X=' num2str(changePoint1) '%'], ...
                            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                    end
                    hold on
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % PD model
                    % from literature for synaptic model
                    AtpRates2_ComplexI = zeros(101,1);
                    cIactivityoptimal2 = zeros(101,1);
                    complexIactivity=model2.ub(ismember(model2.rxns,'NADH2_u10mi'));
                    % change ATP requirement
                    rxn3=model2.rxns(contains(model2.rxns,'ATPM'));
                    model2 = changeRxnBounds(model2,rxn3,ATPMvalue,'l');
                    %                 rxn4=model2.rxns(contains(model2.rxns,'EX_o2[e]'));
                    %                 model2 = changeRxnBounds(model2,rxn4,O2value,'u');
                    % add complex III value???
                    %                 rxn5=model2.rxns(contains(model2.rxns,'CYOR_u10mi'));
                    %                 if any(contains(types{i},'ASYN'))
                    %                     model1 = changeRxnBounds(model1,rxn5,ComplexIII,'u');
                    %                 else
                    %                     model1 = changeRxnBounds(model1,rxn5,ComplexIII,'l');
                    %                 end
                    
                    for n = 0:100;
                        ComplexImodel = model2;
                        % change complex I activity
                        rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
                        ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((n/100)*complexIactivity)),'u');
                        % change objective function
                        rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
                        ComplexImodel = changeObjective(ComplexImodel,rxn2);
                        
                        switch method
                            case 'FBA'
                                FBA_ComplexI = optimizeCbModel(ComplexImodel,'max');
                                FBA_ComplexI.l=FBA_ComplexI.x(find(ismember(ComplexImodel.rxns, rxn1)));
                            case 'eFBA'
                                try
                                    FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);
                                    FBA_ComplexI.f=FBA_ComplexI.v(contains(ComplexImodel.rxns,'ATPS4mi'));
                                    FBA_ComplexI.l=FBA_ComplexI.v(find(ismember(ComplexImodel.rxns, rxn1)));
                                catch
                                    FBA_ComplexI.f=0;
                                    FBA_ComplexI.l=0;
                                end
                        end
                        AtpRates2_ComplexI(n+1) = FBA_ComplexI.f;
                        cIactivityoptimal2(n+1) = FBA_ComplexI.l;
                    end
                    %Plot flux through ATP synthase
                    ratio=AtpRates2_ComplexI/(max(AtpRates2_ComplexI))*100;
                    ratio = round(ratio,2);
                    plot (0:100,ratio, 'LineWidth',2, 'Color',[0 0.8 1]);
                    set(gca, 'FontSize',20);
                    xlabel('Inhibition of Complex I activity (%)');
                    ylabel('Flux through ATP synthase (% of control)');
                    
                    % add anotation
                    hold on;  % Retain the current plot
                    changePoint1=min(find(ratio<99.9)-1);
                    specificY1 = ratio(changePoint1);
                    if ~isempty(changePoint1)
                        plot([changePoint1,changePoint1],[0, specificY1], 'r--', 'LineWidth', 1.5);
                        text(changePoint1,specificY1, ['X=' num2str(changePoint1) '%'], ...
                            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);  % Customize text properties as needed
                    end           %
                    legend(sets{1},'changeline1',sets{2},'changeline2','Box','off','FontSize', 10,'Location', 'best');
                    ylim(gca, [0,108]);
                    xlim(gca, [0,108]);
                    title([sets{1} ' vs ' sets{2} ' with ' method],'FontSize', 14)
                    
                    %Plot flux through Complex I
                    if k==1
                        axes('Position',[.350 .165 .1 .1]);
                    else
                        axes('Position',[.750 .165 .1 .1]);
                    end
                    plot (1:101,cIactivityoptimal1, 'LineWidth',2, 'Color','m');
                    hold on
                    plot (1:101,cIactivityoptimal2, 'LineWidth',2, 'Color',[0 0.8 1]);
                    ylabel(['Complex I flux (NADH2u10mi)';'       (umol/gDW/h)        ']);
                    xlabel('Inhibition of Complex I activity (%)');
                end
            elseif length(modelnum) ==2
                if any(contains(types{i},'ASYN'))
                    ATPMvalue = ATPMvalue_ASYN;
                else
                    ATPMvalue = ATPMvalue_SYN;
                end
                % control and PD models
                model1=All.(types{i}).(modelnum{1});
                model2=All.(types{i}).(modelnum{2});
                
                % control model
                AtpRates1_ComplexI = zeros(11,1);
                cIactivityoptimal1 = zeros(11,1);
                complexIactivity=model1.ub(ismember(model1.rxns,'NADH2_u10mi'));
                % change ATP requirement
                rxn3=model1.rxns(contains(model1.rxns,'ATPM'));
                model1 = changeRxnBounds(model1,rxn3,ATPMvalue,'l');
                num = 0:10:100;
                for m = 1:length(num);
                    ComplexImodel = model1;
                    % change complex I activity
                    rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
                    ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((num(m)/100)*complexIactivity)),'u');
                    % change objective function
                    rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
                    ComplexImodel = changeObjective(ComplexImodel,rxn2);
                    
                    switch method
                        case 'FBA'
                            FBA_ComplexI = optimizeCbModel(ComplexImodel,'max');
                            FBA_ComplexI.l=FBA_ComplexI.x(find(ismember(ComplexImodel.rxns, rxn1)));
                        case 'eFBA'
                            try
                                FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);
                                FBA_ComplexI.f=FBA_ComplexI.v(contains(ComplexImodel.rxns,'ATPS4mi'));
                                FBA_ComplexI.l=FBA_ComplexI.v(find(ismember(ComplexImodel.rxns, rxn1)));
                            catch
                                FBA_ComplexI.f=0;
                                FBA_ComplexI.l=0;
                            end
                    end
                    AtpRates1_ComplexI(m) = FBA_ComplexI.f;
                    cIactivityoptimal1(m) = FBA_ComplexI.l;
                end
                %Plot flux through ATP synthase
                subplot(2,2,i)
                ratio=(AtpRates1_ComplexI/(max(AtpRates1_ComplexI)))*100;
                ratio = round(ratio,2);
                
                % Define a fashionable color palette
                colors = [0, 0.2941, 0.5843;  %  blue
                    0.9255, 0.4784, 0.0314]; % orange
                
                plot (num,ratio, 'LineWidth',2, 'Color',colors(1,:));
                
                % add anotation
                hold on;  % Retain the current plot
                changePoint1=min(find(ratio<99.9)-1);
                specificY1 = ratio(changePoint1);
                specificX1 = num(changePoint1);
                
                hold on
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % PD model
                % from literature for synaptic model
                AtpRates2_ComplexI = zeros(11,1);
                cIactivityoptimal2 = zeros(11,1);
                complexIactivity=model2.ub(ismember(model2.rxns,'NADH2_u10mi'));
                % change ATP requirement
                rxn3=model2.rxns(contains(model2.rxns,'ATPM'));
                model2 = changeRxnBounds(model2,rxn3,ATPMvalue,'l');
                num = 0:10:100;
                for n = 1:length(num);
                    ComplexImodel = model2;
                    % change complex I activity
                    rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
                    ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((num(n)/100)*complexIactivity)),'u');
                    % change objective function
                    rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
                    ComplexImodel = changeObjective(ComplexImodel,rxn2);
                    
                    switch method
                        case 'FBA'
                            FBA_ComplexI = optimizeCbModel(ComplexImodel,'max');
                            FBA_ComplexI.l=FBA_ComplexI.x(find(ismember(ComplexImodel.rxns, rxn1)));
                        case 'eFBA'
                            try
                                FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);
                                FBA_ComplexI.f=FBA_ComplexI.v(contains(ComplexImodel.rxns,'ATPS4mi'));
                                FBA_ComplexI.l=FBA_ComplexI.v(find(ismember(ComplexImodel.rxns, rxn1)));
                            catch
                                FBA_ComplexI.f=0;
                                FBA_ComplexI.l=0;
                            end
                    end
                    AtpRates2_ComplexI(n) = FBA_ComplexI.f;
                    cIactivityoptimal2(n) = FBA_ComplexI.l;
                end
                %Plot flux through ATP synthase
                ratio=AtpRates2_ComplexI/(max(AtpRates2_ComplexI))*100;
                ratio = round(ratio,2);
                plot (num,ratio, 'LineWidth',2, 'Color',colors(2,:));
                set(gca, 'FontSize',10);
                xlabel('Inhibition of Complex I activity (%)','FontSize',14);
                ylabel('Flux through ATP synthase (% of control)','FontSize',14);
                % add anotation for model2
                changePoint2=min(find(ratio<99.9)-1);
                specificY2 = ratio(changePoint2);
                specificX2 = num(changePoint2);
                
                hold on;  % Retain the current plot
                % add anotation for model1
                if ~isempty(changePoint1)
                    plot([specificX1,specificX1],[0, specificY1], 'r--', 'LineWidth', 1.5);
                    text(specificX1,specificY1, ['ChangePoint=' num2str(specificX1) '%'], ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                end
                % add anotation for model2
                if ~isempty(changePoint2)
                    plot([specificX2,specificX2],[0, specificY2], 'r--', 'LineWidth', 1.5);
                    text(specificX2,specificY2, ['ChangePoint=' num2str(specificX2) '%'], ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                end
                ylim(gca, [80,105]);
                yticks(0:50:100);
                xlim(gca, [0,102]);
                title([modelnum{1} ' vs ' modelnum{2} ' with ATPM as ' char(string(ATPMvalue)),' umol/gDW/h'],'FontSize', 15)
                legend(modelnum{1},modelnum{2},'Box','off','FontSize', 12,'Location', 'best');

                %Plot flux through Complex I
                subplot(2,2,i+2)
                %                 if k==1
                %                     axes('Position',[.350 .165 .1 .1]);
                %                 else
                %                     axes('Position',[.750 .165 .1 .1]);
                %                 end
                % Control flux
                plot (num,cIactivityoptimal1, 'LineWidth',2, 'Color',colors(1,:));
                hold on
                % PD flux
                plot (num,cIactivityoptimal2, 'LineWidth',2, 'Color',colors(2,:));
                hold on
                % add changing point
                specificY1 = cIactivityoptimal1(changePoint1);
                if ~isempty(changePoint1)
                    plot([specificX1,specificX1],[0, specificY1], 'r--', 'LineWidth', 1.5);
                    text(specificX1,specificY1, ['ChangePoint=' num2str(specificX1) '%'], ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                end
                specificY2 = cIactivityoptimal2(changePoint2);
                hold on
                if ~isempty(changePoint2)
                    plot([specificX2,specificX2],[0, specificY2], 'r--', 'LineWidth', 1.5);
                    text(specificX2,specificY2, ['ChangePoint=' num2str(specificX2) '%'], ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                end
                
                ylabel(['Complex I flux (NADH2u10mi)';'       (umol/gDW/h)        '],"FontSize",14);
                xlabel('Inhibition of Complex I activity (%)','FontSize',14);
                title(['Complex I flux changes of ' modelnum{1} ' and ' modelnum{2}],'FontSize', 15)
                legend(modelnum{1},modelnum{2},'Box','off','FontSize', 12,'Location', 'best');
            end
        end
    end
end




