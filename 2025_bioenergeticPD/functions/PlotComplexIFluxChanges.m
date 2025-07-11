%% explore the changing trends with different ATPM value
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

types=fieldnames(All);
d = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(types)
    models=All.(types{i});
    modelnum=fieldnames(models);
    % only perform complex I inhibition on unconstrained models
    modelnum=modelnum(~contains(modelnum,'constrain'));% only use unconstrained models!!!
    if length(modelnum) ==2
        % control and PD models
        model1=All.(types{i}).(modelnum{1});
        model2=All.(types{i}).(modelnum{2});
        
        % control model
        AtpRates1_ComplexI = zeros(length(ATPMrange),1);
        cIactivityoptimal1 = zeros(length(ATPMrange),1);
        cIIIactivityoptimal1 = zeros(length(ATPMrange),1);
        for m = 1:length(ATPMrange);
            ATPMvalue = ATPMrange(m);
            ComplexImodel = model1;
            % change ATPM value
            rxn=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
            ComplexImodel = changeRxnBounds(ComplexImodel,rxn,ATPMvalue,'l');
            % change objective function
            rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
            ComplexImodel = changeObjective(ComplexImodel,rxn1);
            % complex I rxn
            rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
            % complex III rxn
            rxn3=ComplexImodel.rxns(contains(ComplexImodel.rxns,'CYOR_u10mi'));
            switch method
                case 'FBA'
                    FBA_ComplexI = optimizeCbModel(ComplexImodel,'max');
                    FBA_ComplexI.l=FBA_ComplexI.x(find(ismember(ComplexImodel.rxns, rxn1)));
                case 'eFBA'
                    try
                        FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);
                        FBA_ComplexI.f=FBA_ComplexI.v(contains(ComplexImodel.rxns,'ATPS4mi'));
                        FBA_ComplexI.l=FBA_ComplexI.v(find(ismember(ComplexImodel.rxns, rxn2)));
                        FBA_ComplexI.c=FBA_ComplexI.v(find(ismember(ComplexImodel.rxns, rxn3)));
                    catch
                        FBA_ComplexI.f=0;
                        FBA_ComplexI.l=0;
                        FBA_ComplexI.c=0;
                    end
            end
            AtpRates1_ComplexI(m) = FBA_ComplexI.f;
            cIactivityoptimal1(m) = FBA_ComplexI.l;
            cIIIactivityoptimal1(m) = FBA_ComplexI.c;
        end
        % PD model
        AtpRates2_ComplexI = zeros(length(ATPMrange),1);
        cIactivityoptimal2 = zeros(length(ATPMrange),1);
        cIIIactivityoptimal2 = zeros(length(ATPMrange),1);
        for n = 1:length(ATPMrange);
            ATPMvalue = ATPMrange(n);
            ComplexImodel = model2;
            % change ATPM value
            rxn=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
            ComplexImodel = changeRxnBounds(ComplexImodel,rxn,ATPMvalue,'l');
            % change objective function
            rxn1=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
            ComplexImodel = changeObjective(ComplexImodel,rxn1);
            % complex I rxn
            rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'NADH2_u10mi'));
            % complex III rxn
            rxn3=ComplexImodel.rxns(contains(ComplexImodel.rxns,'CYOR_u10mi'));
            switch method
                case 'FBA'
                    FBA_ComplexI = optimizeCbModel(ComplexImodel,'max');
                    FBA_ComplexI.l=FBA_ComplexI.x(find(ismember(ComplexImodel.rxns, rxn1)));
                case 'eFBA'
                    try
                        FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);
                        FBA_ComplexI.f=FBA_ComplexI.v(contains(ComplexImodel.rxns,'ATPS4mi'));
                        FBA_ComplexI.l=FBA_ComplexI.v(find(ismember(ComplexImodel.rxns, rxn2)));
                        FBA_ComplexI.c=FBA_ComplexI.v(find(ismember(ComplexImodel.rxns, rxn3)));
                    catch
                        FBA_ComplexI.f=0;
                        FBA_ComplexI.l=0;
                        FBA_ComplexI.c=0;
                    end
            end
            AtpRates2_ComplexI(n) = FBA_ComplexI.f;
            cIactivityoptimal2(n) = FBA_ComplexI.l;
            cIIIactivityoptimal2(n) = FBA_ComplexI.c;
        end
        
        % Define a fashionable color palette
                colors = [0, 0.2941, 0.5843;  % deep blue
                    0.9255, 0.4784, 0.0314; % deep orange
                    0.2196, 0.5059, 0.1843;   % deep green
                    0.5451, 0.7569, 0.9686;   % light blue
                    0.9569, 0.7137, 0.4706;   % light Orange 
                    0.7412, 0.8863, 0.7255];   % light green
        subplot(1,2,i)
        % Control flux
        yyaxis left;
        h1=plot(ATPMrange,cIactivityoptimal1, 'LineWidth', 2, 'Color',colors(1,:));  % Blue line for complexI
        hold on
        h1_1=plot(ATPMrange,cIIIactivityoptimal1, 'LineWidth', 2, 'Color',colors(2,:));
        ylabel('The flux of Complex I and III (umol/gDW/h)');  % Label for the left y-axis
        
        % Plot the second y-axis (right side)
        yyaxis right;
        h2=plot(ATPMrange, AtpRates1_ComplexI, 'LineWidth', 2, 'Color',colors(3,:));  % Red line for ATPS4mi
        ylabel('The flux through ATP synthase (umol/gDW/h)');   % Label for the right y-axis
        hold on
        
        % PD flux
        yyaxis left;
        h3=plot(ATPMrange,cIactivityoptimal2, 'LineWidth', 2, 'Color',colors(4,:));  % Blue line for complexI
        hold on
        h3_1=plot(ATPMrange,cIIIactivityoptimal2, 'LineWidth', 2, 'Color',colors(5,:));
        ylabel('The flux of Complex I and III (umol/gDW/h)','FontSize', 14);  % Label for the left y-axis
 
        
        % Plot the second y-axis (right side)
        yyaxis right;
        h4=plot(ATPMrange, AtpRates2_ComplexI, 'LineWidth', 2, 'Color',colors(6,:));  % Red line for ATPS4mi
        ylabel('The flux through ATP synthase  (umol/gDW/h)','FontSize', 14);   % Label for the right y-axis
        hold on
        legend([h1,h1_1 h2,h3,h3_1,h4], {['Complex I flux in ' modelnum{1}],['Complex III flux in ' modelnum{1}],...
            ['ATP synthase flux in ' modelnum{1}],['Complex I flux in ' modelnum{2}],['Complex III flux in ' modelnum{2}], ['ATP synthase flux in ' modelnum{2}]},'Location', 'best', 'Box','off','FontSize', 12)
        xlabel('Energy demand (ATPM, umol/gDW/h)','FontSize', 14);
        title([modelnum{1} ' and ' modelnum{2}],'FontSize', 15)
    end
end