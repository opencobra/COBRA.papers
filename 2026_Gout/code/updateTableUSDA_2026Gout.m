function fluxTableUsda = updateTableUSDA_2026Gout(nutritionFile, nutritionInfoFile, outputDir)
% Updates the USDA food-flux table with additional/curated nutrition values
% (e.g. purine and urate contents) and saves it as a `.mat` file. The saved
% table is picked up by the Nutrition Toolbox (setInSilicoDiet,
% setFoodRxnsWbm, ...) when the global CBT_NUTRITION_DATABASE_OVERRIDE points
% to outputDir, which this function sets.
%
% USAGE:
%
%    fluxTableUsda = updateTableUSDA_2026Gout(nutritionFile, nutritionInfoFile, outputDir)
%
% INPUTS:
%    nutritionFile:      Path to the additional nutrition values per 100 g
%                        food (e.g. newNutritionValues_gout2026.xlsx). Columns: name,
%                        food ID, database ('usda' or 'frida'), then one
%                        column per VMH metabolite ID.
%    nutritionInfoFile:  Path to the units (ug, mg or g) of every metabolite
%                        column (e.g. newNutritionInfo_gout2026.xlsx), columns metID
%                        and unit.
%
% OPTIONAL INPUTS:
%    outputDir:  Path to the folder where the updated table is written.
%                Defaults to fullfile('updatedDatabases', 'fluxMacroTables').
%
% OUTPUTS:
%    fluxTableUsda:  The updated USDA flux table (mmol per 100 g). Also
%                    written to <outputDir>/USDA<YEAR>_100gFluxValue.mat
%                    (variable `fluxTableUsda`), with <YEAR> the USDA edition
%                    used by the Nutrition Toolbox.
%
% NOTE:
%
%    Nutrition values are converted from ug/mg/g to mmol using the VMH
%    metabolite formulas (getMolecularMass). Only rows with databaseUsed
%    'usda' are used; Frida rows are ignored. A value in the input replaces
%    the database value for that food item and metabolite; empty (NaN) cells
%    leave the database value untouched. Requires the COBRA Toolbox
%    (loadVMHDatabase, getMolecularMass, getNutritionDatabaseFile).
%
% .. Author: - Bram Nap

if nargin < 3 || isempty(outputDir)
    outputDir = fullfile('updatedDatabases', 'fluxMacroTables');
end

% Check the input files
for f = {nutritionFile, nutritionInfoFile}
    if ~isfile(f{1})
        error('updateTableUSDA_2026Gout:missingInput', 'Could not find "%s".', f{1});
    end
end

% Temporarily switch off the database override so the original
% (non-updated) toolbox table is loaded, not a previously saved version.
global CBT_NUTRITION_DATABASE_OVERRIDE
CBT_NUTRITION_DATABASE_OVERRIDE = '';
[baseTablePath, usdaEdition] = getNutritionDatabaseFile('usda', '100gFluxValue');
foodTableUsdaUpdated = load(baseTablePath).fluxTableUsda;

% load new nutrition values
newNutrition = readtable(nutritionFile);
% Load additional information for updated nutrition, such as units
newNutritionInfo = readtable(nutritionInfoFile);

% Extract only USDA items
newNutrition(~strcmpi(newNutrition.databaseUsed, 'usda'),:)= [];

% Metabolite columns start at column 4
metCols = 4:size(newNutrition,2);
metIDs = newNutrition.Properties.VariableNames(metCols);

% get column names = food item ids
cols = foodTableUsdaUpdated.Properties.VariableNames';

% get row names = components
rows = cellstr(foodTableUsdaUpdated.VMHID);

% expand USDA database for new items in components
row = zeros(numel(metIDs),1);
for i = 1:numel(metIDs)
    a = find(strcmpi(rows, metIDs{i}));
    if ~isempty(a)
        row(i) = a;
    else
        foodTableUsdaUpdated.VMHID{end+1,1} = metIDs{i};
        rows{end+1,1} = metIDs{i};
        row(i) = numel(rows);
    end
end

% Convert metabolite values from g to mmol
metabolite2AddValues = newNutrition{:, metCols};

% Load the vmh metabolite database
vmhData = loadVMHDatabase;
vmhMets = vmhData.metabolites;

% Obtain the metabolite formulas. Every metabolite must be found in VMH,
% otherwise the molecular weights would be misaligned with the columns.
[foundMets,~,idxMets] = intersect(metIDs, vmhMets(:,1), 'stable');
if numel(foundMets) ~= numel(metIDs)
    error('updateTableUSDA_2026Gout:unknownMet', 'Metabolite(s) not found in VMH: %s', ...
        strjoin(setdiff(metIDs, foundMets), ', '));
end
met2AddFormulas = vmhMets(idxMets,4);

% Calculate molecular weight
mws = getMolecularMass(met2AddFormulas);

% Find out the order of metabolites of the new nutritional values
[foundInfo,~,metOrder] = intersect(metIDs, newNutritionInfo.metID, 'stable');
if numel(foundInfo) ~= numel(metIDs)
    error('updateTableUSDA_2026Gout:missingUnit', 'No unit given in newNutritionInfo for: %s', ...
        strjoin(setdiff(metIDs, foundInfo), ', '));
end

% Reorder the metabolite info (aka units)
newNutritionInfo = newNutritionInfo(metOrder,:);

% Convert units to conversion factors to grams
nutritionUnit = lower(newNutritionInfo.unit);
nutritionUnit = strrep(nutritionUnit, 'ug', '1e6');
nutritionUnit = strrep(nutritionUnit, 'mg', '1e3');
nutritionUnit = strrep(nutritionUnit, 'g', '1');

% Convert to double
nutritionUnit = str2double(nutritionUnit);
if any(isnan(nutritionUnit))
    error('updateTableUSDA_2026Gout:unknownUnit', 'Units must be ug, mg or g.');
end

% Convert values to g, then mol (g / g/mol) and then mmol
convertedValues = metabolite2AddValues'./nutritionUnit;
convertedValues = convertedValues./mws;
convertedValues = convertedValues'*1000;

% Plug the converted values into the USDA table
for i = 1:size(newNutrition,1)
    databaseItem = newNutrition{i,2};
    if isnumeric(databaseItem)
        databaseItem = string(databaseItem);
    end
    col = find(strcmp(cols, databaseItem));
    if isempty(col)
        warning('updateTableUSDA_2026Gout:foodNotFound', ...
            'USDA food item %s not found in the USDA table and is skipped.', databaseItem);
        continue
    end
    for j = 1:length(row)
        % Only overwrite the database value when a value is given
        if ~isnan(convertedValues(i,j))
            foodTableUsdaUpdated{row(j),col} = convertedValues(i,j);
        end
    end
end

fluxTableUsda = foodTableUsdaUpdated;

% Ensure the output folder exists and save the updated table under the same
% file name as the toolbox table so it can replace it
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
save(fullfile(outputDir, ['USDA', usdaEdition, '_100gFluxValue.mat']), 'fluxTableUsda');

% Point the Nutrition Toolbox to the updated table
CBT_NUTRITION_DATABASE_OVERRIDE = outputDir;
end
