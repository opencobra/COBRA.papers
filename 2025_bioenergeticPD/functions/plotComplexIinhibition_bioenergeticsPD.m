function plotComplexIinhibition_bioenergeticsPD(Allmodels,method)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~isstruct(Allmodels)
    error('check input models')
end

if ~exist('method','var')
    method='FBA';
else
    param.debug = 1;
    param.feasTol=1e-7;
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
    % only perform complex I inhibition on unconstrained models
    modelnum=modelnum(contains(modelnum,'constrain'));% only use unconstrained models!!!
    if ~any(contains(modelnum,'ASTRO'))
        if any(contains(modelnum,'PD'))
            % old vs new
            group.index1=~contains(modelnum,'PD');
            group.index2=contains(modelnum,'PD');
            d = figure('units','normalized','outerposition',[0 0 1 1]);
            if ~any(contains(modelnum,'ASYN'))
                % synaptic info %1.25mg synaptic mitochondrial protein corresponds to 1gWW
                complexIactivity = (136*60*1.25*0.177)/1000;
                oxygenUptake = (58.3*60*1.25*0.177)/1000;
                ATPproduction = (132.6*60*1.25*0.177)/1000;
            else
                %non-synaptic info %5mg non-synaptic mitochondrial protein corresponds to 1gWW
                complexIactivity = (205*60*5*0.177)/1000;
                oxygenUptake = (118*60*5*0.177)/1000;
                ATPproduction = (953*60*5*0.177)/1000;
            end
            for k=1:length(fieldnames(group))
                name=fieldnames(group);
                sets=modelnum(group.(name{k}));
                model1=Allmodels.(types{i}).(sets{1});
                model2=Allmodels.(types{i}).(sets{2});
                if isfield(model1, 'C')
                    model1=rmfield(model1,'C');
                end
                if isfield(model1, 'd')
                    model1=rmfield(model1,'d');
                end
                if isfield(model1, 'ctrs')
                    model1=rmfield(model1,'ctrs');
                end
                if isfield(model1, 'dsense')
                    model1=rmfield(model1,'dsense');
                end
                if isfield(model2, 'C')
                    model2=rmfield(model2,'C');
                end
                if isfield(model2, 'd')
                    model2=rmfield(model2,'d');
                end
                if isfield(model2, 'ctrs')
                    model2=rmfield(model2,'ctrs');
                end
                if isfield(model2, 'dsense')
                    model2=rmfield(model2,'dsense');
                end
                %Generate a model that has all its exchanges closed:
                [selExc,selUpt] = findExcRxns(model1);
                indSyn=find(selUpt);
                uptakeSyn=model1.rxns(indSyn); %uptake reactions
                for n = 1:length(uptakeSyn);
                    model1 = changeRxnBounds(model1,uptakeSyn{n},0,'l');
                end
                % from literature for synaptic model
                AtpRates1_ComplexI = zeros(101,1);
                cIactivityoptimal1 = zeros(101,1);
                for m = 0:100;
                    ComplexImodel = model1;
                    if ~any(contains(modelnum,'ASYN'))
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',-10,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',10,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',-10,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',10,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_adn[e]',-100,'l');
%                         if strcmp(method , 'eFBA');
%                             ComplexImodel = changeRxnBounds(ComplexImodel,'EX_cholate[e]',-9.69606,'l');
%                             ComplexImodel = changeRxnBounds(ComplexImodel,'EX_cholate[e]',9.69606,'u');
%                         end
                    else
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',-1,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',1,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',-1,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',1,'u');
                    end
                    ComplexImodel = changeRxnBounds(ComplexImodel,'EX_o2[e]',-oxygenUptake,'l');
                    ComplexImodel = changeRxnBounds(ComplexImodel,'ATPM',ATPproduction, 'u');
                    rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
%                     if strcmp(method , 'eFBA')
%                         ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((m/100)*complexIactivity)),'b');
%                     else
                        ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((m/100)*complexIactivity)),'u');
%                     end
                  if ~isfield(ComplexImodel,'c') | strcmp(method,'FBA') | (isfield(ComplexImodel,'c') & ~any(ComplexImodel.c))
                        rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
                        % rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
                        ComplexImodel = changeObjective(ComplexImodel,rxn2);
                      end
                    switch method
                        case 'FBA'
                            FBA_ComplexI = optimizeCbModel(ComplexImodel,'max');
                            if FBA_ComplexI.stat ==1
                            FBA_ComplexI.l=FBA_ComplexI.x(find(ismember(ComplexImodel.rxns, rxn1)));
                            else
                                FBA_ComplexI.l=0;
                                FBA_ComplexI.f=0;
                            end
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
                if ~any(AtpRates1_ComplexI)
                    error('model infeasible')
                end
                %Plot flux through ATP synthase
                subplot(1,2,k)
                ratio=(AtpRates1_ComplexI/(max(AtpRates1_ComplexI)))*100;
                plot (1:101,ratio, 'LineWidth',2, 'Color',[1 0 1]);
                
                % add anotation
                hold on;  % Retain the current plot
                changePoint=min(find(ratio<99.99)-1);
                specificY = ratio(changePoint);
                %                 plot(changePoint,specificY,'ro', 'MarkerSize', 10);
                plot([changePoint,changePoint],[0, specificY], 'r--', 'LineWidth', 1.5);
                text(changePoint,specificY, ['X=' num2str(changePoint) '%'], ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                hold on
                %
                % add another plots
                [selExc,selUpt] = findExcRxns(model2);
                indSyn=find(selUpt);
                uptakeSyn=model2.rxns(indSyn); %uptake reactions
                for n = 1:length(uptakeSyn);
                    model2 = changeRxnBounds(model2,uptakeSyn{n},0,'l');
                end
                % from literature for synaptic model
                AtpRates2_ComplexI = zeros(101,1);
                cIactivityoptimal2 = zeros(101,1);
                for n = 0:100
                    ComplexImodel = model2;
                    if ~any(contains(modelnum,'ASYN'))
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',-10,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',10,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',-10,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',10,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_adn[e]',-100,'l');
%                         if strcmp(method , 'eFBA');
%                             ComplexImodel = changeRxnBounds(ComplexImodel,'EX_cholate[e]',-9.69606,'l');
%                             ComplexImodel = changeRxnBounds(ComplexImodel,'EX_cholate[e]',9.69606,'u');
%                         end
                    else
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',-1,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',1,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',-1,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',1,'u');
                    end
                    ComplexImodel = changeRxnBounds(ComplexImodel,'EX_o2[e]',-oxygenUptake,'l');
                    ComplexImodel = changeRxnBounds(ComplexImodel,'ATPM',ATPproduction, 'u');
                    rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
%                     if strcmp(method , 'eFBA')
%                         ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((n/100)*complexIactivity)),'b');
%                     else
                        ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((n/100)*complexIactivity)),'u');
%                     end
                      if ~isfield(ComplexImodel,'c') | strcmp(method,'FBA') | (isfield(ComplexImodel,'c') & ~any(ComplexImodel.c))
                        rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
                        % rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
                        ComplexImodel = changeObjective(ComplexImodel,rxn2);
                      end
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
                plot (1:101,ratio, 'LineWidth',2, 'Color',[0 0.8 1]);
                set(gca, 'FontSize',20);
                xlabel('Inhibition of Complex I activity (%)');
                ylabel('Flux through ATP synthase (% of control)');
                
                % add anotation
                hold on;  % Retain the current plot
                changePoint=min(find(ratio<99.9)-1);
                specificY = ratio(changePoint);
                plot([changePoint,changePoint],[0, specificY], 'r--', 'LineWidth', 1.5,'Color','blue');
                text(changePoint,specificY, ['X=' num2str(changePoint) '%'], ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                %
                legend(sets{1},'changeline1',sets{2},'changeline2','Box','off','FontSize', 10);
                ylim(gca, [0,108]);
                xlim(gca, [0,108]);
                title([sets{1} ' vs ' sets{2} ' with ' method],'FontSize', 15)
                
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %      PD vs control
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            group.index1=contains(modelnum,'1');
            group.index2=contains(modelnum,'2');
            d = figure('units','normalized','outerposition',[0 0 1 1]);
            if ~any(contains(modelnum,'ASYN'))
                % synaptic info
                complexIactivity = (136*60*1.25*0.177)/1000;
                oxygenUptake = (58.3*60*1.25*0.177)/1000;
                ATPproduction = (132.6*60*1.25*0.177)/1000;
            else
                %non-synaptic info
                complexIactivity = (205*60*5*0.177)/1000;
                oxygenUptake = (118*60*5*0.177)/1000;
                ATPproduction = (953*60*5*0.177)/1000;
            end
            for k=1:length(fieldnames(group))
                name=fieldnames(group);
                sets=modelnum(group.(name{k}));
                model1=Allmodels.(types{i}).(sets{1});
                model2=Allmodels.(types{i}).(sets{2});
                if isfield(model1, 'C')
                    model1=rmfield(model1,'C');
                end
                if isfield(model1, 'd')
                    model1=rmfield(model1,'d');
                end
                if isfield(model1, 'ctrs')
                    model1=rmfield(model1,'ctrs');
                end
                if isfield(model1, 'dsence')
                    model1=rmfield(model1,'dsence');
                end
                if isfield(model2, 'C')
                    model2=rmfield(model2,'C');
                end
                if isfield(model2, 'd')
                    model2=rmfield(model2,'d');
                end
                if isfield(model2, 'ctrs')
                    model2=rmfield(model2,'ctrs');
                end
                if isfield(model2, 'dsence')
                    model2=rmfield(model2,'dsence');
                end
                %Generate a model that has all its exchanges closed:
                [selExc,selUpt] = findExcRxns(model1);
                indSyn=find(selUpt);
                uptakeSyn=model1.rxns(indSyn); %uptake reactions
                for n = 1:length(uptakeSyn);
                    model1 = changeRxnBounds(model1,uptakeSyn{n},0,'l');
                end
                % from literature for synaptic model
                AtpRates1_ComplexI = zeros(101,1);
                cIactivityoptimal1 = zeros(101,1);
                for m = 0:100;
                    ComplexImodel = model1;
                    if ~any(contains(modelnum,'ASYN'))
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',-10,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',10,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',-10,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',10,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_adn[e]',-100,'l');
%                         if strcmp(method , 'eFBA');
%                             ComplexImodel = changeRxnBounds(ComplexImodel,'EX_cholate[e]',-9.69606,'l');
%                             ComplexImodel = changeRxnBounds(ComplexImodel,'EX_cholate[e]',9.69606,'u');
%                         end
                    else
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',-1,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',1,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',-1,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',1,'u');
                    end
                    ComplexImodel = changeRxnBounds(ComplexImodel,'EX_o2[e]',-oxygenUptake,'l');
                    ComplexImodel = changeRxnBounds(ComplexImodel,'ATPM',ATPproduction, 'u');
                    rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
%                     if strcmp(method , 'eFBA')
%                         ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((m/100)*complexIactivity)),'b');
%                     else
                        ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((m/100)*complexIactivity)),'u');
%                     end
                      if ~isfield(ComplexImodel,'c') | strcmp(method,'FBA') | (isfield(ComplexImodel,'c') & ~any(ComplexImodel.c))
                        rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
                        % rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
                        ComplexImodel = changeObjective(ComplexImodel,rxn2);
                      end
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
                plot (1:101,ratio, 'LineWidth',2, 'Color',[1 0 1]);
                
                % add anotation
                hold on;  % Retain the current plot
                changePoint=min(find(ratio<99.99)-1);
                specificY = ratio(changePoint);
                %                 plot(changePoint,specificY,'ro', 'MarkerSize', 10);
                plot([changePoint,changePoint],[0, specificY], 'r--', 'LineWidth', 1.5);
                text(changePoint,specificY, ['X=' num2str(changePoint) '%'], ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                hold on
                %
                % add another plots
                [selExc,selUpt] = findExcRxns(model2);
                indSyn=find(selUpt);
                uptakeSyn=model2.rxns(indSyn); %uptake reactions
                for n = 1:length(uptakeSyn);
                    model2 = changeRxnBounds(model2,uptakeSyn{n},0,'l');
                end
                % from literature for synaptic model
                AtpRates2_ComplexI = zeros(101,1);
                cIactivityoptimal2 = zeros(101,1);
                for n = 0:100;
                    ComplexImodel = model2;
                    if ~any(contains(modelnum,'ASYN'))
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',-10,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',10,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',-10,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',10,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_adn[e]',-100,'l');
%                         if strcmp(method , 'eFBA');
%                             ComplexImodel = changeRxnBounds(ComplexImodel,'EX_cholate[e]',-9.69606,'l');
%                             ComplexImodel = changeRxnBounds(ComplexImodel,'EX_cholate[e]',9.69606,'u');
%                         end
                    else
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',-1,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_glu_L[e]',1,'u');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',-1,'l');
                        ComplexImodel = changeRxnBounds(ComplexImodel,'EX_mal_L[e]',1,'u');
                    end
                    ComplexImodel = changeRxnBounds(ComplexImodel,'EX_o2[e]',-oxygenUptake,'l');
                    ComplexImodel = changeRxnBounds(ComplexImodel,'ATPM',ATPproduction, 'u');
                    rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
%                     if strcmp(method , 'eFBA')
%                         ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((n/100)*complexIactivity)),'b');
%                     else
                        ComplexImodel = changeRxnBounds(ComplexImodel,rxn1,(complexIactivity - ((n/100)*complexIactivity)),'u');
%                     end
                      if ~isfield(ComplexImodel,'c') | strcmp(method,'FBA') | (isfield(ComplexImodel,'c') & ~any(ComplexImodel.c))
                        rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
                        % rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
                        ComplexImodel = changeObjective(ComplexImodel,rxn2);
                      end
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
                plot (0:100,ratio, 'LineWidth',2, 'Color',[0 0.8 1]);
                set(gca, 'FontSize',20);
                xlabel('Inhibition of Complex I activity (%)');
                ylabel('Flux through ATP synthase (% of control)');
                
                % add anotation
                hold on;  % Retain the current plot
                changePoint=min(find(ratio<99.9)-1);
                specificY = ratio(changePoint);
                plot([changePoint,changePoint],[0, specificY], 'r--', 'LineWidth', 1.5,'Color','blue');
                text(changePoint,specificY, ['X=' num2str(changePoint) '%'], ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);  % Customize text properties as needed
                %
                legend(sets{1},'changeline1',sets{2},'changeline2','Box','off','FontSize', 10);
                ylim(gca, [0,108]);
                xlim(gca, [0,108]);
                title([sets{1} ' vs ' sets{2} ' with ' method],'FontSize', 15)
                
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
            
        end
    end
    
end

