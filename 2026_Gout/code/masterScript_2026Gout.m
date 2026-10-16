% Master script for performed for the paper:
% Urate metabolism in gout: predicting the effects of diet and the gut 
% microbiome using whole-body metabolic modelling

% Builds food-item whole-body models
% (germ-free Harvey/Harvetta and microbiome WBMs), predicts the maximal blood
% urate flux on pre-defined diets and on individual food items, identifies
% the microbial contributors to these fluxes and plots the results.
%
% Run this script from the folder it is in (2026_Gout/code). The input
% files are read from 2026_Gout/inputData. All results are written below
% resultDir. The functions in this folder are put on top of the MATLAB path
% so they take precedence over older copies.
%
% Requires the COBRA Toolbox (with the Nutrition Modelling Toolbox and the
% PSCM toolbox), an LP solver, the Statistics and Machine Learning Toolbox
% (lasso, bootci) and the Parallel Computing Toolbox.
%
% .. Author: - Bram Nap

%% User input
% Folder with the microbiome WBMs (mWBMs, .mat files). Model file names must
% contain _<ID>_ with <ID> as listed in inputData/metadataForDietSetting.xlsx.
mWBMDir = '';
% Relative abundance table of the microbes in each mWBM (rows are mWBM IDs,
% columns are taxa), used to find the microbial contributors and clusters
% in steps 7 and 8
microbeAbundancePath = '';

% Folder where all results are saved. A relative path is relative to this
% folder.
resultDir = 'results';

% Number of parallel workers
numWorkers = 8;

%% Setup: check input and set all paths
codeDir = pwd;
if ~isfile(fullfile(codeDir, 'driver_dietTestingUrate.m'))
    error('Please run this driver from the folder it is located in (2026_Gout/code).')
end
if isempty(mWBMDir) || ~isfolder(mWBMDir)
    error('Please set mWBMDir to the folder with the mWBM .mat files.')
end
if isempty(microbeAbundancePath) || ~isfile(microbeAbundancePath)
    error('Please set microbeAbundancePath to the relative abundance table of the mWBMs.')
end

% Ensure the functions in this folder take precedence
addpath(codeDir);

% Make resultDir absolute so functions that use what() on a folder cannot
% pick up a different folder with the same name on the MATLAB path
if isempty(regexp(resultDir, '^([A-Za-z]:|[\\/])', 'once'))
    resultDir = fullfile(codeDir, resultDir);
end

% Input files shipped in the inputData folder (next to the code folder)
inputDir = fullfile(fileparts(codeDir), 'inputData');
diet2Set = fullfile(inputDir, 'dietFormulationForDietSetting.xlsx'); % Diet template: one food item, replaced in the FBAs
metadataFile = fullfile(inputDir, 'metadataForDietSetting.xlsx'); % Links each mWBM ID to its sex and the diet template
foodMetadata = fullfile(inputDir, 'shortNameGoutFoodsForFigures.xlsx'); % Short names and Eat/Limit/Avoid category of the food items
% Additional purine/urate food composition values per 100 g food and the
% units of each metabolite column
nutritionValuesFile = fullfile(inputDir, 'newNutritionValues_gout2026.xlsx');
nutritionInfoFile = fullfile(inputDir, 'newNutritionInfo_gout2026.xlsx');

% Check that all input files are present
for inputFile = {diet2Set, metadataFile, foodMetadata, nutritionValuesFile, nutritionInfoFile}
    if ~isfile(inputFile{1})
        error('Input file not found: %s', inputFile{1});
    end
end

% All result paths branch off resultDir
updatedDatabaseDir = fullfile(resultDir, 'updatedDatabase'); % Supplemented USDA table
gfFoodDir = fullfile(resultDir, 'foodItemWBMs', 'germFree'); % Germ-free food-item WBMs (step 2)
mwbmFoodDir = fullfile(resultDir, 'foodItemWBMs', 'mWBM'); % Microbiome food-item WBMs (step 3)
fbaDir = fullfile(resultDir, 'fba'); % FBA results (steps 4 and 5)
microbeDir = fullfile(resultDir, 'microbialContributions'); % Microbial contributors and clusters (steps 7 and 8)
figureDir = fullfile(resultDir, 'figures'); % Figures (step 9)

% setInSilicoDiet saves the food-item WBMs in <outputDir>/foodDietWBMs.
% These folders are used as model input for the FBAs.
modelDirs = struct('germFree', fullfile(gfFoodDir, 'foodDietWBMs'), ...
    'mWBM', fullfile(mwbmFoodDir, 'foodDietWBMs'));

% Metabolite to maximise
metabolite = {'urate[bc]'};

% FBA runs: which food-item WBMs are used, whether pre-defined diets or
% individual food items are tested and whether diets are normalised to
% 2000 kcal. Results of each run go to <fbaDir>/<name>.
fbaRuns = struct( ...
    'name',      {'existingDietGF', 'existingDietGFNorm', 'existingDietMWBM', 'foodItemGF', 'foodItemMWBM'}, ...
    'models',    {'germFree',       'germFree',           'mWBM',             'germFree',   'mWBM'}, ...
    'test',      {'diet',           'diet',               'diet',             'food',       'food'}, ...
    'normalise', {false,            true,                 false,              false,        false});

%% Step 1 Supplement the USDA database with the additional nutrition values
% The updated table is saved in updatedDatabaseDir and the global
% CBT_NUTRITION_DATABASE_OVERRIDE is set to it, so the Nutrition Toolbox
% uses it instead of the original USDA table in steps 2 and 3.
updateTableUSDA_2026Gout(nutritionValuesFile, nutritionInfoFile, updatedDatabaseDir);

%% Step 2 Set food items on Harvey and Harvetta (germ-free)
% Ensure the supplemented USDA table is used (also when this section is run
% on its own)
global CBT_NUTRITION_DATABASE_OVERRIDE
CBT_NUTRITION_DATABASE_OVERRIDE = updatedDatabaseDir;

setInSilicoDiet(diet2Set, 'outputDir', gfFoodDir, 'constrainFluxWBM', false, ...
    'constrainFoodWBM', true, 'checkFeasibility', false);

% Allow dietary urate to be absorbed from the gut into the blood (by default
% it only passes to the feces). The food-item WBMs are overwritten.
addUrateUptake2WBMs_2026Gout(modelDirs.germFree, numWorkers);

%% Step 3 Set food items on the mWBMs
global CBT_NUTRITION_DATABASE_OVERRIDE
CBT_NUTRITION_DATABASE_OVERRIDE = updatedDatabaseDir;

setInSilicoDiet(diet2Set, 'outputDir', mwbmFoodDir, 'constrainFluxWBM', false, ...
    'constrainFoodWBM', true, 'checkFeasibility', false, 'metadataPath', metadataFile, ...
    'pathToWbms', mWBMDir);

% Allow dietary urate to be absorbed from the gut into the blood. The
% food-item WBMs are overwritten.
addUrateUptake2WBMs_2026Gout(modelDirs.mWBM, numWorkers);

%% Step 4 Test pre-existing diets for urate production
% NOTE: the FBAs use every .mat file in the food-item WBM folders of steps 2
% and 3. Remove models from earlier runs from these folders first.
existingDiet2Test = {
    'DACHDiet.txt'
    'EUAverageDietNew.m'
    'GlutenFreeDiet.txt'
    'HighFatLowCarbDiet.txt'
    'HighFiberDiet.m'
    'HighProteinDiet.m'
    'MediterraneanDiet.txt'
    'Type2DiabetesDiet.txt'
    'UnhealthyDiet.m'
    'VeganDiet.txt'
    'VegetarianDiet.m'
    };

% Caloric value (kcal) of each diet, in the same order as existingDiet2Test.
% Used to normalise the diets to 2000 kcal.
dietKcal = [2148.36
    2371.9
    2983.83
    2215.48
    2435.67
    1791.84
    2647.45
    2084.57
    3847.05
    1991.21
    1888.29];

saveFullRes = true;

for r = find(strcmp({fbaRuns.test}, 'diet'))
    runDir = fullfile(fbaDir, fbaRuns(r).name);
    if ~isfolder(runDir); mkdir(runDir); end
    testFluxExistingDiet_2026Gout(existingDiet2Test, metabolite, runDir, modelDirs.(fbaRuns(r).models), ...
        saveFullRes, fbaRuns(r).normalise, dietKcal);
end

%% Step 5 Test individual contributions of food items
% USDA food item IDs, each tested at 100 g
foodItems = {'173589'
    '173430'
    '171267'
    '168568'
    '170457'
    '330137'
    '170112'
    '172183'
    '171287'
    '173032'
    '169118'
    '167762'
    '173944'
    '169998'
    '170272'
    '173800'
    '169291'
    '170124'
    '169708'
    '170844'
    '171706'
    '172818'
    '172449'
    '169736'
    '174284'
    '168389'
    '173741'
    '169986'
    '172473'
    '167723'
    '170567'
    '169251'
    '167914'
    '172451'
    '170427'
    '170379'
    '170826'
    '167895'
    '171966'
    '172850'
    '175132'
    '174037'
    '172934'
    '172584'
    '174210'
    '171077'
    '172922'
    '175167'
    '174727'
    '172420'
    '174729'
    '173709'
    '173717'
    '174182'};

saveFullRes = true;

for r = find(strcmp({fbaRuns.test}, 'food'))
    runDir = fullfile(fbaDir, fbaRuns(r).name);
    if ~isfolder(runDir); mkdir(runDir); end
    testFluxFoodItems_2026Gout(foodItems, metabolite, runDir, modelDirs.(fbaRuns(r).models), ...
        saveFullRes, '', numWorkers);
end

%% Step 6 Test individual dietary metabolites for urate production
% (runCalcIndMet) For every food-item WBM: maximal blood urate flux with
% 100 mmol of each individual dietary metabolite added on top of the EU
% average diet, using calcIndMetContrFromFood_2026Gout.

% All dietary metabolites in the (supplemented) USDA and the Frida table
global CBT_NUTRITION_DATABASE_OVERRIDE
CBT_NUTRITION_DATABASE_OVERRIDE = updatedDatabaseDir;
load(getNutritionDatabaseFile('usda', '100gFluxValue'), 'fluxTableUsda');
load(getNutritionDatabaseFile('frida', '100gFluxValue'), 'fluxTableFrida');

dietMetList = unique([fluxTableUsda.VMHID; fluxTableFrida.VMHID]);
% Add the purines and urate only if they are not already in the list
dietMetList = [dietMetList; setdiff({'urate';'ade'; 'gua';'xan'; 'hxan'}, dietMetList, 'stable')];

indMetDir = fullfile(resultDir, 'indDietMets');
for modelSet = {'germFree', 'mWBM'}
    modelPath = modelDirs.(modelSet{1});
    model = what(modelPath).mat;
    for i = 1:size(model,1)
        modRes = calcIndMetContrFromFood_2026Gout(fullfile(modelPath, model{i}), dietMetList, metabolite{1});

        if i == 1
            totalRes = modRes;
        else
            totalRes(:,i+1) = modRes(:,2);
        end
    end

    % Header row: MetID and the model names (used by the figure in step 9)
    totalRes = [[{'MetID'}, erase(model', '.mat')]; totalRes];

    saveDir = fullfile(indMetDir, modelSet{1});
    if ~isfolder(saveDir); mkdir(saveDir); end
    writecell(totalRes, fullfile(saveDir, 'dietMetaboliteContributions.xlsx'));
    clear totalRes
end

%% Step 7 Find the microbial contributors to the flux values (mWBM runs)
% Assumes a single metabolite: one flux table per run is analysed.
if ~isscalar(metabolite)
    error('Step 7 is set up for a single metabolite.')
end

% Parameters for identifying the most important microbial species
param = struct;
param.bootSamp = 10000; % Boot samples for obtaining the mean and 95%CI of microbe contribution potentials
param.minFreq = 0.5; % Minimal lasso selection frequency
param.nBootLasso = 500; % Boot samples for the lasso stability selection
param.enBoot = 100000; % Boot samples for the elastic-net regressions

% Start parallel pool if needed
poolobj = gcp('nocreate'); if isempty(poolobj); parpool(numWorkers); end

% Diet labels as used in the result tables (file extension removed)
dietLabels = regexprep(existingDiet2Test, '\.(m|txt)$', '');

for r = find(strcmp({fbaRuns.models}, 'mWBM'))
    runDir = fullfile(fbaDir, fbaRuns(r).name);
    contrDir = fullfile(microbeDir, fbaRuns(r).name);
    if ~isfolder(contrDir); mkdir(contrDir); end

    % Reaction names (<diet or food item>_<metabolite>) and the flux table
    if strcmp(fbaRuns(r).test, 'diet')
        param.rxnsOfInterest = strcat(dietLabels, '_', metabolite{1});
        fluxFilePath = fullfile(runDir, ['fbaResDiet_', metabolite{1}, '.xlsx']);
    else
        param.rxnsOfInterest = strcat(foodItems, '_', metabolite{1});
        fluxFilePath = fullfile(runDir, ['fbaResFood_', metabolite{1}, '.xlsx']);
    end

    % Extract the microbe biomass shadow prices from the FBA solutions
    extractMicrobeContributions_2026Gout(runDir, contrDir);

    % Find the associated microbial contributors to the fluxes
    findMicrobialContributors_2026Gout(fullfile(contrDir, 'biomass_shadow_prices'), fluxFilePath, ...
        microbeAbundancePath, contrDir, param);
end

%% Step 8 Find the microbial clusters that correlate best with the fluxes
% For each flux, clusters of 1-5 of its microbial contributors (step 7) are
% correlated (Spearman) with the flux, using the summed relative abundance
% of the cluster. The best positive and negative cluster per cluster size is
% saved in <microbeDir>/<run>/clusterCorrelations.xlsx.
for r = find(strcmp({fbaRuns.models}, 'mWBM'))
    runDir = fullfile(fbaDir, fbaRuns(r).name);
    contrDir = fullfile(microbeDir, fbaRuns(r).name);
    if strcmp(fbaRuns(r).test, 'diet')
        fluxFilePath = fullfile(runDir, ['fbaResDiet_', metabolite{1}, '.xlsx']);
    else
        fluxFilePath = fullfile(runDir, ['fbaResFood_', metabolite{1}, '.xlsx']);
    end
    findBestMicrobeClusters_2026Gout(fluxFilePath, microbeAbundancePath, ...
        fullfile(contrDir, 'enetFluxMicrobeRes.csv'), contrDir);
end

%% Step 9 Plot the results
% Germ-free runs (one female and one male model) are shown as bar charts,
% mWBM runs as box plots. Food item runs use the food short names and the
% Eat/Limit/Avoid categories. The individual metabolite results of step 6
% are shown as heatmaps.
if ~isfolder(figureDir); mkdir(figureDir); end
ylabelText = {'Predicted urate flux in circulating blood'; '(mmol/person/day)'};

for r = 1:numel(fbaRuns)
    runDir = fullfile(fbaDir, fbaRuns(r).name);
    if strcmp(fbaRuns(r).test, 'diet')
        dataFilePath = fullfile(runDir, ['fbaResDiet_', metabolite{1}, '.xlsx']);
        metadataForPlot = '';
    else
        dataFilePath = fullfile(runDir, ['fbaResFood_', metabolite{1}, '.xlsx']);
        metadataForPlot = foodMetadata;
    end

    if strcmp(fbaRuns(r).models, 'germFree')
        barchartNutritionSimulations_2026Gout(dataFilePath, figureDir, fbaRuns(r).name, ylabelText, ...
            metadataForPlot, '', metabolite{1});
    else
        boxplotNutritionSimulation_2026Gout(dataFilePath, figureDir, fbaRuns(r).name, ylabelText, ...
            metadataForPlot, '', metabolite{1});
    end
end

% Increase in urate production per individual dietary metabolite (step 6):
% germ-free and average mWBM per sex, and per mWBM for selected metabolites
indMetDir = fullfile(resultDir, 'indDietMets');
plotIndMetContribution_2026Gout(fullfile(indMetDir, 'germFree', 'dietMetaboliteContributions.xlsx'), ...
    fullfile(indMetDir, 'mWBM', 'dietMetaboliteContributions.xlsx'), figureDir, 'indMetContribution', ...
    {'5fthf', 'cholp'});
