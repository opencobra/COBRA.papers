% Sampling model & uncertainty reduction
% Uniform sampling of feasible solution space, was used to predict a set of 
% metabolites that could potentially be secreted by the iDopaNeuro models. 
% In principle, all of these metabolites should be considered for targeted 
% metabolomic analyses in exometabolomic experiments.

% RHMC or CHRR is called via the function sampleCbModel. The main inputs to
% sampleCbModel are a COBRA model structure, the name of the selected
% sampler and a parameter struct that controls properties of the sampler
% used. In the instance of CHRR, two parameters are important: the sampling
% density (nStepsPerPoint) and the number of samples (nPointsReturned). The
% total length of the random walk is nStepsPerPoint*nPointsReturned. The
% time it takes to run the sampler depends on the total length of the
% random walk and the size of the model. k

% For CHRR, an additional on/off parameter (toRound) controls whether or not the
% polytope is rounded. Rounding large models can be slow but is strongly
% recommended for the first round of sampling. Below we show how to get
% around this step in subsequent rounds. The method outputs two results.
% First, the rounded sample space of the model used for sampling (in case
% of toRound = 1 this would be the rounded model), and second, the samples
% generated.

clear

% Select a COBRA model

pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

% Path to additional images for Uniform sampling figure (experimental design)
%i.e. covariance.png  eDistance.png  fluxSpace.png
auxImagesDir = [pathSave filesep 'auxImages'];
    
% Select the sampler
% optionsSampling.samplerName = 'CHRR';
optionsSampling.samplerName = 'RHMC';
    
% Select solver
[~, ~] = changeCobraSolver('mosek', 'all', 0);

models = {'iDopaNeuroCT'; 'iDopaNeuroC'};
for k = 1:length(models)
    
    % Load model
    load([pathSave filesep models{k} filesep models{k} '.mat'])
    eval(['model = ' models{k} ';'])
    
    % Sample the flux space
    if isfile([pathSave filesep models{k} filesep models{k} '_samples_' optionsSampling.samplerName '.mat'])
        fprintf('%s\n',[models{k} '_samples_' optionsSampling.samplerName '.mat found, loaded.'])
        load([pathSave filesep models{k} filesep models{k} '_samples_' optionsSampling.samplerName '.mat'])
        optionsSampling.nPointsReturned = size(samples,2);
    else
        fprintf('%s\n',[models{k} '_samples_' optionsSampling.samplerName '.mat not found, being computed.'])
        switch optionsSampling.samplerName
            case 'CHRR'
                
                % Round the polytope to calculate the nPointsReturned
                if ~isfile([pathSave filesep models{k} filesep 'P_model.mat'])
                    
                    % The polytopes are rounded and the sampling points are ingnored.
                    % Rounding flux space
                    optionsSampling.nStepsPerPoint = 1;
                    optionsSampling.nPointsReturned = 1;
                    optionsSampling.toRound = 1;
                    
                    [P_model, ~] =  sampleCbModel(model, [], 'CHRR', optionsSampling);
                    save([pathSave filesep  models{k} filesep 'P_model.mat'], 'P_model')
                    
                else
                    load([pathSave filesep models{k} filesep 'P_model.mat'])
                end
                optionsSampling.nPointsReturned = 8 * size(P_model.N, 2);
                
                % Depending the size of the polytopes  the optimal sampling its is
                % selected as  and . This time, we can avoid the rounding step by inputting
                % the rounded polytope from the previous round of sampling (this may take
                % some time).
                tic
                optionsSampling.toRound = 0;
                optionsSampling.nStepsPerPoint = 1000; %size(P_model.N, 2)^2;
                %[modelSampling,samples,volume] = sampleCbModel(model, sampleFile, samplerName, options, modelSampling)
                [~, samples] =  sampleCbModel(model, [], optionsSampling.samplerName, optionsSampling, P_model);
                optionsSampling.nPointsReturned = size(samples,2);
                toc
                
            case 'RHMC'
                tic
                %[modelSampling,samples,volume] = sampleCbModel(model, sampleFile, samplerName, options, modelSampling)
                [~, samples] =  sampleCbModel(model, [], optionsSampling.samplerName, optionsSampling);
                toc
        end
        save([pathSave filesep models{k} filesep models{k} '_samples_' optionsSampling.samplerName '.mat'], 'samples')
    end
    
    
    %% Calculate covariance matrix and obtain TopKNorms
    covS = cov(samples');
    
    % Exchange metabolite table
    if 0
        nRows = numOfExchanges;
        varTypes = {'double', 'string', 'string', 'double', 'double', 'double', 'double'};
        varNames = {'rank', 'rxns', 'rxnNames', 'lb', 'ub', 'S_norms', 'S_adjusted_norms'};
        exMetTable = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes,...
            'VariableNames', varNames);
        rxnIdxs = findRxnIDs(model, model.rxns(S_rxns));
        exMetTable.rank = (1:numOfExchanges)';
        exMetTable.rxns = model.rxns(rxnIdxs);
        exMetTable.rxnNames = model.rxnNames(rxnIdxs);
        exMetTable.lb = model.lb(rxnIdxs);
        exMetTable.ub = model.ub(rxnIdxs);
        exMetTable.S_norms = S_adjusted_norms;
        exMetTable.S_adjusted_norms = S_norms;
    end
    
    % Set the number of top reactions
    %     numOfExchanges = sum(~model.SIntRxnBool) - 20;
    numOfExchanges = 20;
    
    %[A_rxns, adjustedA_norms, A_norms] = getTopKNormsAdaptive(A,k,SIntRxnBool)
    [S_rxns,S_adjusted_norms, S_norms] = getTopKNormsAdaptive(covS, numOfExchanges, model.SIntRxnBool);
    plot_data = [S_adjusted_norms S_norms - S_adjusted_norms];
    
    %% Figure
    % Obtain the names of the metabolites with the greatest degree of uncertainty.
    met_names = {};
    for j = 1:numOfExchanges
        met_names{j} = model.metNames{find(full(model.S(:,strcmp(model.rxns{S_rxns(j)},model.rxns))))};
        mets{j} = model.mets{find(full(model.S(:,strcmp(model.rxns{S_rxns(j)},model.rxns))))};
    end
    [~, ia, ~] = unique(met_names);
    met_names = met_names(sort(ia));
    plot_data = plot_data(sort(ia), :);
    logPlot_data = [log10(plot_data(:, 1))  log10(sum(plot_data, 2)) - log10(plot_data(:, 1))];
    met_names = regexprep(met_names, 'Smcfa-Blood-Pool', 'Short/medium chain fatty acids');
    met_names = regexprep(met_names, 'All Trans ', '');
       
    %%
    param.objectives = {'unWeightedTCBMflux'};
    solutions = modelMultipleObjectives(model, param);
    
    %% Uniform sampling figure (experimental design)
    %figure('units','normalized','outerposition',[0 0 1 1]) % one screen
    figure('units','normalized','outerposition',[0 0 0.5 1]) %two screens
    % Sampling
    subplot(12, 5, 5 * (1:4) - 4);
    imshow([auxImagesDir filesep 'fluxSpace.png'])
    title('A. Uniform sampling', 'FontSize', 14)
    % Covariance matrix
    subplot(12, 5, 5 * (5:8) - 4);
    imshow([auxImagesDir filesep 'covariance.png'])
    title('B. Covariance matrix', 'FontSize', 14)
    % Euclidean norms
    subplot(12, 5, 5 * (9:12) - 4);
    imshow([auxImagesDir filesep 'eDistance.png'])
    title('C. Euclidean norms', 'FontSize', 14)
    % Priority metabolites
    subplot(12, 5, 5 * (1:12) - 2);
    barh(logPlot_data,'stacked');
    xlabel('$log_{10} (Euclidean \: Norm)$', 'interpreter', 'latex')
    title('D. Priority metabolites', 'FontSize', 14)
    set(gca,'yticklabels',S_rxns);
    yticks(1:length(met_names))
    set(gca, 'yticklabel', regexprep(met_names,'\_', '\\_'), 'YDir','reverse', 'FontSize', 14);
    % ax = gca;
    % % This sets background color to magenta
    % ax.YColor = 'magenta';
    
    % Sampled distributions
    rxnsToPlot = model.rxns(S_rxns);
    for i = 1:4
        switch i
            case 1
                slots = [4 9 14 19 24];
                idxFigure = 'E';
            case 2
                slots = slots + 1;
                idxFigure = 'F';
            case 3
                slots = [39 44 49 54 59];
                idxFigure = 'G';
            case 4
                slots = slots + 1;
                idxFigure = 'H';
        end
        subplot(12, 5, slots);
        currRxn = findRxnIDs(model, rxnsToPlot{i});
        [f, xi] = hist(samples(currRxn, :), 100);
        %H = histogram(samples(currRxn, :), 100);
        
        f = f / sum(f); % frequency to probablility
        % Flux distribution
        area(xi, f,'FaceColor', 'r','FaceAlpha', .2, 'EdgeAlpha', .2);
        hold on
        % Mean
        scatter(mean(samples(currRxn, :)), max(f), 'filled');
        
        % Standard deviation
        errorbar(mean(samples(currRxn, :)),max(f),std(samples(currRxn, :)), ...
            'horizontal');
        
        % Prediction
        if solutions.unWeightedTCBMflux.stat == 1
            xline(solutions.unWeightedTCBMflux.v(currRxn), '-b', 'LineWidth', 2);
        end
        % legend('p1', 'sdD', 'meanD', 'Min flux FVA', 'Max flux FVA', 'prediction', 'location', 'best')
        title({[idxFigure '. Flux distribution of'], met_names{i}}, 'FontSize', 14)
        xlabel({'Flux distribution', '(µmol/gDW/hr)'}, 'FontSize', 14)
        ylabel('Probability', 'FontSize', 14)
        
        % axis
        axis([min([xi solutions.unWeightedTCBMflux.v(currRxn)]) - ((max(xi) - min(xi)) * 0.1) ...
            max([xi solutions.unWeightedTCBMflux.v(currRxn)]) + ((max(xi) - min(xi)) * 0.1) ...
            0 max(f) + ((max(f) - min(f)) * 0.1)])
        
        if i == 1
            legend('Marginal distribution', 'Sample Mean', 'Sample standard deviation', ...
                'entropy maximisation (flux)', 'location', 'southoutside', 'FontSize', 12)
        end
    end

    
    %% save the figure
    
    
    savefig([pathSave filesep models{k} filesep 'priorityMetabolites_' optionsSampling.samplerName ...
        '_' int2str(optionsSampling.nPointsReturned) '_' models{k}])
    saveas(gcf,[pathSave filesep models{k} filesep 'priorityMetabolites_' optionsSampling.samplerName ...
        '_' int2str(optionsSampling.nPointsReturned) '_' models{k} '.eps']);
    saveas(gcf,[pathSave filesep models{k} filesep 'priorityMetabolites_' optionsSampling.samplerName ...
        '_' int2str(optionsSampling.nPointsReturned) '_' models{k} '.png']);
end