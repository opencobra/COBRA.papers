function [elasticNetResults, elasticNetStats, lassoResTab, microbeContributionStats] = findMicrobialContributors_2026Gout(mContributionDir, fluxPath, mappedMicrobePath, saveDir, param)
% Identifies the microbial taxa most associated with a set of flux predictions
% by combining bootstrapped microbe biomass shadow-price statistics with lasso
% stability selection and bootstrapped elastic-net regression. Cleaned up
% version of microbeToFluxPipeline with only the parts necessary for its
% functioning.
%
% USAGE:
%
%    [elasticNetResults, elasticNetStats, lassoResTab, microbeContributionStats] = findMicrobialContributors_2026Gout(mContributionDir, fluxPath, mappedMicrobePath, saveDir, param)
%
% INPUTS:
%    mContributionDir:      Path to the directory with the microbial biomass
%                           shadow prices associated with the flux predictions.
%    fluxPath:              Path to the file with processed flux predictions.
%    mappedMicrobePath:     Path to the relative abundance table of all
%                           microbial taxa in each model.
%    saveDir:               Path to the directory where all result tables are
%                           saved.
%    param:                 Structure with fields:
%
%                             * .rxnsOfInterest - list of reactions to investigate
%                             * .bootSamp - number of bootstrap samples for the 95% CI around the shadow-price means
%                             * .nBootLasso - number of bootstrap samples for the lasso stability selection
%                             * .minFreq - minimal lasso stability-selection frequency for further analysis
%                             * .enBoot - number of bootstrap samples for the elastic-net regressions (used to generate per-coefficient statistics: 95% CIs, standard errors, t-values, p-values)
%
% OUTPUTS:
%    elasticNetResults:         Elastic-net regression coefficients of the
%                               microbial predictors for the associated flux
%                               predictions.
%    elasticNetStats:           Summary statistics of the elastic-net
%                               regression model fit.
%    lassoResTab:               Selection frequency of each predictor variable
%                               (e.g. microbial taxa) in the bootstrapped lasso
%                               regressions, i.e. the fraction of bootstrap
%                               samples for which a non-zero coefficient was
%                               obtained.
%    microbeContributionStats:  Table with the mean and bootstrapped 95%
%                               confidence intervals of the flux-associated
%                               shadow prices of the microbial taxa in the
%                               model.
%
% .. Author: - Tim Hensen, July 2025 - April 2026

% Set random number generator for reproducible results
rng(1, "twister")

% Extract the function parameters
rxnsOfInterest = param.rxnsOfInterest;
bootSamp = param.bootSamp;
enBoot = param.enBoot;
minFreq = param.minFreq;
nBootLasso = param.nBootLasso;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Identifying potential flux influencers in the microbiomes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the microbial species with the highest average flux contribution
% potentials
poolobj = gcp('nocreate');
numRxns = numel(rxnsOfInterest);
if isempty(poolobj); parpool(numRxns); end

% Perform bootstrap analysis on all rxnsOfInterest
bootMeanTables = cell(1,numRxns);
disp('Shadow price analysis')
parfor i = 1:numRxns
    bootMeanTables{i} = getMicrobeSensitivity(mContributionDir, rxnsOfInterest{i}, bootSamp);
end

% Create tall table with 
bootMeanTable = vertcat(bootMeanTables{:});

% Remove microbial species with an average shadow price of zero within the
% 95% CI.
bootMeanTableFiltered = bootMeanTable( bootMeanTable.('p>0.05')==0 , :);

% Generate intersection table 

% Get microbe-metabolite pairs
metMicrobePairs = bootMeanTableFiltered(:,{'Metabolite','Taxa'}); 

% Create intersection table with taxa names as row names
metMicrobePairs.value = ones(height(metMicrobePairs),1);
metMicrobePairsWide = unstack(metMicrobePairs,"value","Metabolite",'VariableNamingRule','preserve');
metMicrobePairsWide.Properties.RowNames = metMicrobePairsWide.Taxa;

% Replace nan with zero
fluxLimTable = fillmissing(metMicrobePairsWide(:,2:end),'constant',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LASSO feature selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Filter the fluxes table on the reactions in fluxLimTable
% fluxLimTable = readtable(fluxLimTablePath,'VariableNamingRule','preserve','ReadRowNames',true);
rxnsOfInterest = fluxLimTable.Properties.VariableNames;
vmhIDs = renameVmhToMetName_2026Gout(rxnsOfInterest, true);

% Load the processed fluxes and process sample IDs
fluxes = readtable(fluxPath,'VariableNamingRule','preserve');
% The fbaRes*.xlsx tables written by testFluxExistingDiet_2026Gout/testFluxFoodItems_2026Gout
% store the model IDs in a first column called 'Row'
if ~any(strcmp(fluxes.Properties.VariableNames, 'ID')) && any(strcmp(fluxes.Properties.VariableNames, 'Row'))
    fluxes = renamevars(fluxes, 'Row', 'ID');
end
% Remove quote characters around IDs and sex (present in manually edited tables)
fluxes.ID = erase(fluxes.ID, {'''', '"'});
fluxes.Sex = erase(fluxes.Sex, {'''', '"'});
fluxes = fluxes(:,[{'ID','Sex'}, vmhIDs]); % Filter on metabolites of interest
fluxes.Properties.VariableNames = [{'ID','Sex'}, rxnsOfInterest];
fluxes.ID = strrep(fluxes.ID, 'foodItem', '');
fluxes.ID = strrep(fluxes.ID, '.mat', '');


% Load the relative abundances
microbiome = readtable(mappedMicrobePath,'VariableNamingRule','preserve','ReadRowNames',false);
microbiome = renamevars(microbiome,'Row','ID');    
if isvar(microbiome,'Sum of taxa')
    microbiome = removevars(microbiome,'Sum of taxa');
end

% Set all nan values to zero
microbiome(:,2:end) = fillmissing(microbiome(:,2:end),'constant',0);

% Ensure that the samples are in the correct order between the fluxes and
% relative abundances
[~,ia,ib] = intersect(fluxes.ID,microbiome.ID,'stable');
fluxes = fluxes(ia,:);
microbiome = microbiome(ib,:);

% Get the associated microbial species for each reaction
fluxMicrobes = cellfun(@(x) fluxLimTable.Row(fluxLimTable.(x)==1),rxnsOfInterest,'UniformOutput',false);

% Combine tables into a single table with the fluxes at the end
mergeFun = @(x,y) outerjoin(microbiome(:,[{'ID'},y']),fluxes(:,{'ID','Sex',x}),'Keys','ID','MergeKeys',true,'Type','left');
fluxMicrobeData = cellfun(mergeFun,rxnsOfInterest,fluxMicrobes,'UniformOutput',false);

% Get the number of microbial species with any contribution to the fluxes
microbeContributionStats = cell(1,2);
microbeContributionStats{1} = getMeanSDNumMicrobes(fluxMicrobeData);

% Remove the ID column
fluxMicrobeData = cellfun(@(x) removevars(x,'ID'), fluxMicrobeData,'UniformOutput',false); 

% Find subset of microbial species that predict the associated fluxes
disp('Feature selection')
lassoRes = pruneMicrobialFeaturesWithLasso(fluxMicrobeData, nBootLasso);

% Create table for output
lassoResTab = vertcat(lassoRes{:});

% Create wide table for the lasso selection frequencies
lassoResTab = unstack(lassoResTab,"selection_frequency","Reaction",'VariableNamingRule','preserve');

% Filter microbes with a frequency below minFreq
lassoResFiltered = cellfun(@(x) x(x.selection_frequency>minFreq,:), lassoRes,'UniformOutput',false);

% Prune microbes from relative abundance data
microbeNames = cellfun(@(x) x.Properties.VariableNames(1:end-2), fluxMicrobeData,'UniformOutput',false); % Get the microbe names
microbesToKeep = cellfun(@(x) x.Species,lassoResFiltered,'UniformOutput',false);

% Remove microbes 
pruneMicrobes = @(x,y,z) removevars(x,setdiff(y,z) );
fluxMicrobeData = cellfun(pruneMicrobes,fluxMicrobeData, microbeNames ,microbesToKeep,'UniformOutput',false);

% Get the number of microbial contributors again after lasso-based pruning
microbeContributionStats{2} = getMeanSDNumMicrobes(fluxMicrobeData);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Feature importance calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify relevant features
[enetTab, elasticNetStats] = elasticNetMicrobeWeights(fluxMicrobeData, enBoot);

elasticNetResults = vertcat(enetTab{:});
elasticNetResults = renamevars(elasticNetResults,'Metabolite','Reaction');
elasticNetResults = addvars(elasticNetResults, renameVmhToMetName_2026Gout(elasticNetResults.Reaction),'NewVariableNames','Metabolite','After','Reaction');

% Convert elastic net results to wide table with coefficients for each
% flux-microbe association
% enetTabWide = unstack(vertcat(enetTab{:}),'Beta','Metabolite','VariableNamingRule','preserve');
% enetTabWide(:,2:end) = fillmissing(enetTabWide(:,2:end),'constant',0);
% elasticNetResults = enetTabWide;

% Convert VMH IDs to metabolite names
% elasticNetResults.Properties.VariableNames(2:end) = renameVmhToMetName_2026Gout(elasticNetResults.Properties.VariableNames(2:end));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save the shadow price averages for each flux-microbe pair
filePath = fullfile(saveDir,'bootciMicrobeSPvalues.csv');
writetable(bootMeanTable,filePath,'WriteRowNames',true)

% Save intersection table
fluxLimTablePath = fullfile(saveDir,'fluxMicrobeAssociations.csv');
writetable(fluxLimTable,fluxLimTablePath,'WriteMode','overwrite','WriteRowNames',true)

% Save lasso estimation results
filePath = fullfile(saveDir,'lassoMicrobeSelectionFreq.csv');
writetable(lassoResTab,filePath,'WriteRowNames',true)

% Save elastic net estimation results
filePath = fullfile(saveDir,'enetFluxMicrobeRes.csv');
writetable(elasticNetResults,filePath,'WriteRowNames',true)

% Save summary statistics
filePath = fullfile(saveDir,'fluxMicrobeSummaryStats.xlsx');
writetable(microbeContributionStats{1},filePath,'WriteRowNames',true,'Sheet','Contributors');
writetable(microbeContributionStats{2},filePath,'WriteRowNames',true,'Sheet','Analysed');
writetable(elasticNetStats,filePath,'WriteRowNames',true,'Sheet','EnetStats');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bootMeanTable = getMicrobeSensitivity(mContributionDir, metabolite, bootSamp)


% Load flux contributions (files are saved as <rxn>.csv)
contributionFile = fullfile(mContributionDir, metabolite);
if ~isfile(contributionFile)
    contributionFile = [contributionFile, '.csv'];
end
contributions = readtable(contributionFile, 'PreserveVariableNames', true,'ReadRowNames',true);

rmPan = @(x) renamevars(x,x.Properties.VariableNames,erase(x.Properties.VariableNames,'pan'));
contributions = rmPan(contributions);

% Set all nan values to zeros
contributions = fillmissing(contributions, 'constant', 0);

% Filter on microbial species that are associated with the fluxes 
% in one or more samples
contributions = contributions(:,any(abs(contributions{:,:})>1e-6));

% Find the average flux contributions 
contributionMatrix = table2array(contributions);

% Find microbial species without any potential contributions
nzIdx = any(contributionMatrix~=0);

% Reduce contribution matrix
contributionMatrix = contributionMatrix(:,nzIdx);

% Calculate the mean flux contributions and the associated 95% confidence
% intervals
[ci,bootstat] = arrayfun(@(x) bootci(bootSamp,{@mean,contributionMatrix(:,x)},'Alpha',0.05), 1:size(contributions,2),'UniformOutput',false);

% Get means and CI
bootMeans = [cellfun(@mean,bootstat)',[ci{:}]'];

% Check if the 95% CI crossses zero
decimal = 4;
bootMeans(:,4) = ( round(bootMeans(:,2),decimal) .* round(bootMeans(:,3),decimal) ) <= 0;

% Generate table
bootMeanTable = array2table(bootMeans,'VariableNames',{'Mean','2.5CI','97.5CI','p>0.05'});

% Add taxa
bootMeanTable.Taxa = contributions.Properties.VariableNames(nzIdx)';

% Add associated metabolite
bootMeanTable.Metabolite = cellstr(repmat(metabolite,height(bootMeanTable),1));

% Process table
bootMeanTable = movevars(bootMeanTable,'Metabolite','before',1);
bootMeanTable = movevars(bootMeanTable,'Taxa','after','Metabolite');
end


function lassoRes = pruneMicrobialFeaturesWithLasso(fluxMicrobeData, nBootLasso)
%
% INPUTS:
% fluxMicrobeData       Cell array with the microbe relative 
%                       abundances from the WBMs, sample sex 
%                       information, and a final column with 
%                       the associated flux predictions.
%
% OUTPUTS:
% lassoResFiltered      Cell array with a list of microbial 
%                       species for which non-zero beta coefficients 
%                       were obtained when applying a lambda 
%                       regularization estimator that resulted 
%                       in the smallest mean square error (MSE) 
%                       for the flux prediction, while controlling 
%                       for sex.


% Encode sex information into binary variable (0/1) and apply function to
% each cell.
cf = @(f,x) cellfun(f, x, 'UniformOutput',false); % Preallocate cellfun function for compact code
convertSex = @(x) convertvars(x,'Sex', @(x) grp2idx(categorical(x))-1);
microbeDataCell = cf(convertSex, fluxMicrobeData);

% Get predictors and response variables and z-scaled the data column-wise
transforDat = @(x) normalize(table2array(x));
predictors = cf(@(x) transforDat(x(:,1:end-1)), microbeDataCell);
responses = cf(@(x) transforDat(x(:,end)), microbeDataCell);

% Fetch microbe and metabolite names
microbeNames = cf(@(x) x.Properties.VariableNames(1:end-1), fluxMicrobeData); % Get the microbe names
metNames = cf(@(x) x.Properties.VariableNames(end), fluxMicrobeData); % Get the metabolite names

% Preallocate replicated metabolite names for lasso table annotation
metNameForTable = cellfun(@(x,y) repmat(x,length(y),1), metNames, microbeNames,'UniformOutput',false);

% Preallocate table for annotated selection frequencies
lassoRes = cellfun(@(x,y) table(x, y','VariableNames',{'Reaction','Species'}),metNameForTable,microbeNames, 'UniformOutput',false);

% Perform lasso regression on each cell
numMets = numel(responses);
numSamp = height(fluxMicrobeData{1});
tic
fprintf("> Microbe selection in progress... \n")
for i = 1:numMets

    % Get predictor and response vars for metabolite i
    pred = predictors{i};
    resp = responses{i};

    % Create selection of random samples that take 70% of samples each time
    rsampIdx = arrayfun(@(x) randsample(numSamp, round(0.7 * numSamp),true), 1:nBootLasso, 'UniformOutput', false);
    
    % Generate samples for predictors
    X_sub = cellfun(@(x) pred(x, :), rsampIdx,'UniformOutput', false);
    y_sub = cellfun(@(x) resp(x, :), rsampIdx,'UniformOutput', false);
    
    % Preallocate result vector
    selection_count = zeros(size(pred,2), 1); 
    parfor j = 1:nBootLasso
        % Perform lasso regression on random sample
        [b,lStat] = lasso(X_sub{j},y_sub{j},'CV',10,'Intercept',false,'Standardize',false); % Data is already normalized
    
        % Use the lambda with minimum CV error
        coef = b(:, lStat.IndexMinMSE);
        
        % Record which variables were selected
        selection_count = selection_count + (coef ~= 0);
    end
    
    % Add the selection frequencies to table
    lassoRes{i}.selection_frequency = selection_count / nBootLasso;

    % Display progress
    disp(append("> ",string(round(i/numMets*100,0)),"%"))
end
toc

end

function meanSdNumMicrobesForMets = getMeanSDNumMicrobes(fluxMicrobeData)
% This function obtains the mean average and SD of microbial contributors.
% 
% INPUTS:
% fluxMicrobeData           Cell array with tables with gut microbiome data
%
% OUTPUTS:
% meanSdNumMicrobesForMets  Table with the mean average number of microbes
%                           present for each metabolite.

if matches('ID',fluxMicrobeData{1}.Properties.VariableNames)
    fluxMicrobeData = cellfun(@(x) removevars(x,'ID'), fluxMicrobeData,'UniformOutput',false); 
end
microbeNames = cellfun(@(x) x.Properties.VariableNames(1:end-2), fluxMicrobeData,'UniformOutput',false);  % Get microbe names
metNames = cellfun(@(x) x.Properties.VariableNames(end), fluxMicrobeData,'UniformOutput',true); % Get the metabolite names

% Filter on microbes in the microbiome
getMicrobesFun = @(x,y) table2array(x(:,y'));
subsets = cellfun(getMicrobesFun,fluxMicrobeData, microbeNames,'UniformOutput',false);

% Get the number of microbes per sample
msdFun = @(x) [sum(any(x,1)), mean( sum(x~=0,2) ), std( sum(x~=0,2) )];
meanSdForMets = cellfun(msdFun,subsets,'UniformOutput',false);
meanSdForMets = vertcat(meanSdForMets{:});
meanSdNumMicrobesForMets = array2table(meanSdForMets,'VariableNames',{'Total','Mean','SD'},'RowNames',metNames');
end


function [enetTab, enet_stats] = elasticNetMicrobeWeights(fluxMicrobeData, enBoot)
% This function finds the weight coefficients of microbial species for flux
% predictions using elastic net regression

alpha = 0.5; % Elastic net ridge parameters

% Encode sex information into binary variable (0/1) and apply function to
% each cell.
cf = @(f,x) cellfun(f, x, 'UniformOutput',false); % Preallocate cellfun function for compact code
convertSex = @(x) convertvars(x,'Sex', @(x) grp2idx(categorical(x))-1);
microbeDataCell = cf(convertSex, fluxMicrobeData);

% Get predictors and response variables and z-scaled the data column-wise
transforDat = @(x) normalize(table2array(x));
predictors = cf(@(x) transforDat(x(:,1:end-1)), microbeDataCell);
responses = cf(@(x) transforDat(x(:,end)), microbeDataCell);

% Fetch microbe and metabolite names
microbeNames = cf(@(x) x.Properties.VariableNames(1:end-1), fluxMicrobeData); % Get the microbe names
metNames = cf(@(x) x.Properties.VariableNames(end), fluxMicrobeData); % Get the metabolite names

% Preallocate replicated metabolite names for lasso table annotation
metNameForTable = cellfun(@(x,y) repmat(x,length(y),1), metNames, microbeNames,'UniformOutput',false);

% Preallocate table for elastic net results
enetTab = cellfun(@(x,y) table(x, y', zeros(length(y'),1),'VariableNames',{'Metabolite','Taxa','Beta'}),metNameForTable, microbeNames, 'UniformOutput',false);

% Preallocate table with summary statistics
numMets = length(metNames);
varNames = {'MSE','RMSE', 'SST', 'SSR','R2','DoF','Alpha','minLambda'};
% enet_stats = array2table(zeros(numMets,length(varNames)),'VariableNames',varNames,'RowNames',string(metNames)');

enet_stats = cellfun(@(x) array2table(zeros(1,length(varNames)),'VariableNames',varNames,'RowNames',x), metNames, 'UniformOutput',false);

% Warm up algorithm
lasso(predictors{1}, responses{1}, 'CV',10,'Intercept',false,'Standardize',false,'Alpha',alpha);
lasso(predictors{2}, responses{2}, 'CV',10,'Intercept',false,'Standardize',false,'Alpha',alpha);

parfor i=1:numMets
    % Get predictor and response vars for metabolite i
    pred = predictors{i};
    resp = responses{i};
    
    % Remove nan values
    nanSamples = isnan(resp);
    y = resp(~nanSamples);
    X = pred(~nanSamples,:);

    if 1 
        [coef, se, ci, t_stat, p_val, MSE, DF, LambdaMinMSE, SST, SSR, R2] = perfBootEnet(X, y , enBoot);

        % Add regression statistics
        enetTab{i}.Beta = coef;
        enetTab{i} = addvars(enetTab{i}, ci(:,1), ci(:,2), se, t_stat, p_val,  repmat(enBoot, length(coef),1), 'NewVariableNames', {'2.5% CI','97.5% CI', 'SE','tStat','pVal','enBoot'});

        enet_stats{i}{1,:} = [MSE, sqrt(MSE), SST, SSR, R2, DF, alpha, LambdaMinMSE];
        % enet_stats.DoF(i) = DF; % Degrees of freedom
        % enet_stats.minLambda(i) = LambdaMinMSE;
        % enet_stats.R2(i) = R2;
    else 
        % Perform lasso regression on random sample
        % [b,lStat] = lasso(X, y, 'CV',10,'Intercept',false,'Standardize',false,'Alpha',alpha); % Data is already normalized
        % 
        % % Find the coefficients for each microbe and add them to the table
        % coef = b(:, lStat.IndexMinMSE);
        % enetTab{i}.Beta = coef;
        % 
        % % Populate statistics table
        % enet_stats.MSE(i) = lStat.MSE(lStat.IndexMinMSE); % MSE
        % enet_stats.DoF(i) = lStat.DF(lStat.IndexMinMSE); % Degrees of freedom
        % enet_stats.minLambda(i) = lStat.LambdaMinMSE;
        % 
        % % Add the R2 value
        % y_pred = X * coef + lStat.Intercept(lStat.IndexMinMSE); % Predict using the selected coefficients and intercept
        % SST = sum((y - mean(y)).^2);         % total sum of squares (SST) in y
        % SSR = sum((y - y_pred).^2);          % residual sum of squares (SSR) variance
        % enet_stats.R2(i) = 1 - SSR / SST;   % R2
    end
end

enet_stats = vertcat(enet_stats{:});

% Populate alpha variable
% enet_stats.Alpha = repmat(alpha,height(enet_stats),1);

end

function [coef, se, ci, t_stat, p_val, MSE, DF, LambdaMinMSE, SST, SSR, R2] = perfBootEnet(X, y , enBoot)
% Perform bootstrapped elastic net regressions


% Preallocate matrices
n = length(y);
bootCoefs = zeros(size(X,2), enBoot);
regStats = zeros(enBoot,4);

% Generate bootstrap samples
rsampIdx = arrayfun(@(x) randsample(n, round(1*n,0),true), 1:enBoot, 'UniformOutput', false);
X_sub = cellfun(@(x) X(x, :), rsampIdx,'UniformOutput', false);
y_sub = cellfun(@(x) y(x, :), rsampIdx,'UniformOutput', false);

fprintf('Bootstrapping...');
for i = 1:enBoot
    % Run Elastic Net on bootstrap sample
    [b, stat] = lasso(X_sub{i}, y_sub{i}, 'CV',10,'Intercept',false,'Standardize',false,'Alpha',0.5); 
    bootCoefs(:,i) = b(:, stat.IndexMinMSE);
    regStats(i,:) = [stat.Intercept(stat.IndexMinMSE), stat.DF(stat.IndexMinMSE), stat.MSE(stat.IndexMinMSE), stat.LambdaMinMSE];
end
fprintf('Done.\n');

% --- Calculate Statistics ---
DF = mean(regStats(:,2));
MSE = mean(regStats(:,3));
LambdaMinMSE = mean(regStats(:,4));

% Calculate the R2 value
coef = mean(bootCoefs,2);
y_pred = X * coef;                   % Predict using the selected coefficients and intercept
SST = sum((y - mean(y)).^2);         % total sum of squares (SST) in y
SSR = sum((y - y_pred).^2);          % residual sum of squares (SSR) variance
R2 = 1 - SSR / SST;                  % R2

se = std(bootCoefs, 0, 2);                      % Standard Error
ci = prctile(bootCoefs, [2.5 97.5], 2);         % 95% CI Lower
t_stat = coef ./ se;                            % T-statistic (approximate)

% P-value: Proportion of times the bootstrap coef is 0 (or crosses zero)
p_val = 2 * min(mean(bootCoefs > 0, 2), mean(bootCoefs < 0, 2));
end
