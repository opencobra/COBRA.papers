function testFluxFoodItems_2026Gout(foodItems, metabolite, resultDir, modelPath, saveFullRes, foodItemsFile, numWorkers)
% Maximises the production of one or more metabolites for each individual food
% item across every whole-body metabolic (WBM) model in a directory, and writes
% the objective values (and optionally the full FBA solutions) to disk.
%
% USAGE:
%
%    testFluxFoodItems_2026Gout(foodItems, metabolite, resultDir, modelPath, saveFullRes)
%    testFluxFoodItems_2026Gout(foodItems, metabolite, resultDir, modelPath, saveFullRes, foodItemsFile)
%    testFluxFoodItems_2026Gout(foodItems, metabolite, resultDir, modelPath, saveFullRes, foodItemsFile, numWorkers)
%
% INPUTS:
%    foodItems:      Cell array of food item identifiers (USDA IDs) matching
%                    the `Food_EX_<id>_usda` exchange reactions in the models.
%                    Each item is added at 100 g uptake from the USDA database.
%                    Ignored when a valid `foodItemsFile` is supplied.
%    metabolite:     Cell array of metabolite identifiers including the
%                    compartment (e.g. `'urate[bc]'`) to maximise through a
%                    demand reaction.
%    resultDir:      Char, directory where the result files are written.
%    modelPath:      Char, directory containing the WBM models (`.mat`).
%    saveFullRes:    Logical, if true the full FBA solution (fluxes `v`,
%                    shadow prices `y`, reduced costs `w` and microbe biomass
%                    shadow prices) is saved for every model.
%
% OPTIONAL INPUTS:
%    foodItemsFile:  Char, path to a headerless, three-column text file that
%                    defines the food items to test, one per row:
%                       column 1 - food item ID number
%                       column 2 - source database: 'usda' or 'frida'
%                                  (selects the `_usda`/`_frida` reaction suffix)
%                       column 3 - uptake amount in grams
%                    When given and the file exists it overrides `foodItems`
%                    (and the per-item database and amount). When omitted, empty
%                    or pointing at a missing file, reading is skipped and the
%                    `foodItems` argument is used (USDA ids at 100 g).
%    numWorkers:     Numeric, number of workers used by the `parfor` loop over
%                    the WBM models. When > 0 a parallel pool of this size is
%                    started (an existing pool of a different size is
%                    recreated); when 0 the loop runs serially on the client
%                    (useful for debugging). Defaults to 1.
%
% OUTPUTS:
%    None. Files are written to `resultDir`:
%
%                     * fbaSol_<model>.mat - full FBA solution per model (when saveFullRes is true)
%                     * fbaResFood_<metabolite>.xlsx - objective value per model (row) and food item (column)
%
% NOTE:
%
%    Each food item is tested in isolation: all `Food_EX_` reactions are first
%    closed, then the single item is opened to its uptake amount. Dietary
%    urate uptake (gut to blood) is not added here; it is added to the
%    food-item WBMs once, by addUrateUptake2WBMs_2026Gout in the driver.
%    Requires the COBRA Toolbox and the WBM/Persephone helper functions
%    (addDemandReaction, changeRxnBounds, changeObjective, optimizeWBModel).
%
% .. Author: - Bram Nap

% Default the optional food-items file to empty (skip reading) when not given.
if nargin < 6
    foodItemsFile = '';
end

% Default the number of parallel workers to 1 when not given.
if nargin < 7 || isempty(numWorkers)
    numWorkers = 1;
end

% Per-item database suffix and uptake amount. The defaults reproduce the
% original behaviour: every item is a USDA item taken up at 100 g.
foodItems  = foodItems(:);
foodDb     = repmat({'usda'}, numel(foodItems), 1);
foodAmount = repmat(100, numel(foodItems), 1);

% If a food-items file is defined AND present, read the per-item id (column 1),
% source database (column 2 -> _usda/_frida suffix) and uptake amount in grams
% (column 3) from it, overriding the foodItems argument. If the file is not
% defined or not found, reading is skipped and foodItems is used as-is.
if ~isempty(foodItemsFile) && isfile(foodItemsFile)
    foodTable  = readtable(foodItemsFile, 'ReadVariableNames', false);
    foodItems  = cellstr(string(foodTable{:, 1}));    % ids as char vectors
    foodDb     = lower(cellstr(string(foodTable{:, 2})));
    foodAmount = foodTable{:, 3};
elseif ~isempty(foodItemsFile)
    warning('testFluxFoodItems_2026Gout:foodFileMissing', ...
        'Food-items file not found (%s); using the foodItems argument instead.', foodItemsFile);
end

% Find the file names of the wbm models to be tested
wbmModels = what(modelPath).mat;

% Initialise the result variable
results = zeros([length(wbmModels), length(foodItems), length(metabolite)]);
% Initialise the sex variable (preallocated so it can be sliced by parfor).
% Kept as a 1-by-N row to match the original grow-from-{} orientation, since
% the table assignment below transposes it (sex').
sex = cell(1, length(wbmModels));

% Define rxn IDs, one per solution column (food item x metabolite)
if isscalar(foodItems)
    % If only 1 food item is tested the metabolites are unique ids
    fullRxnList = metabolite(:);
else
    % If multiple food items are tested, add food item id to metabolite
    % name for unique rxn ids (<food>_<metabolite>, same as the result table
    % columns).
    fullRxnList = cell([length(metabolite)*length(foodItems),1]);
    idxCount = 0;
    for k = 1:length(foodItems)
        % allocate the names in the list
        startIdx = idxCount*length(metabolite) + 1;
        endIdx = startIdx + length(metabolite) - 1;
        fullRxnList(startIdx:endIdx,1) = strcat(foodItems(k), '_', metabolite);
        % Increase counter
        idxCount = idxCount + 1;
    end
    % Remove .mat and .txt from the names
    fullRxnList = strrep(fullRxnList, '.m', '');
    fullRxnList = strrep(fullRxnList, '.txt', '');
end

% Capture the currently configured LP solver so each parallel worker can
% re-select it (workers start with a clean COBRA environment). Following the
% same pattern as analyseWBMs.m, the parent-process environment (paths and
% globals) is captured once and restored on every worker.
global CBT_LP_SOLVER
solver = CBT_LP_SOLVER;
environment = getEnvironment();

% Set up the parallel pool with the requested number of workers. An existing
% pool of the same size is reused; a differently sized pool is recreated. When
% numWorkers is 0 the pool is left untouched and the parfor below runs serially
% on the client. Mirrors the pool handling in analyseWBMs.m.
if numWorkers > 0
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers);
    elseif poolobj.NumWorkers ~= numWorkers
        delete(poolobj);
        parpool(numWorkers);
    end
end

parfor (i = 1:length(wbmModels), numWorkers)
    % Restore the COBRA environment and solver on this worker
    restoreEnvironment(environment);
    if ~isempty(solver)
        changeCobraSolver(solver, 'LP', 0, -1);
    end

    % Load the model
    model2Test = load(fullfile(modelPath, wbmModels{i}));
    % Store the sex of the model
    sex{i} = model2Test.sex;

    model2Test = setDietConstraints(model2Test, 'EUAverageDiet');

    % Add demand reactions for metabolites that currently do not have it
    for z = 1:size(metabolite,2)
        if ~any(strcmp(model2Test.rxns, strcat('DM_', metabolite(z))))          
            model2Test = addDemandReaction(model2Test, metabolite(z));
        end
    end
    % Show which model number is currently being investigated
    disp(i)

    % Defined in every iteration (filled only when saveFullRes is true) so
    % parfor does not treat it as an uninitialised temporary
    solutionStruct = struct();

    if saveFullRes
        % Initialise solution variable if the FBA struct has to be saved
        solutionStruct = struct();
        
        % Set the rxn ids to the solution structure
        solutionStruct.rxns = fullRxnList;
        
        % Add model names
        solutionStruct.ID = erase(wbmModels(i),'.mat');
        
        % Add model sex
        solutionStruct.sex = cellstr(model2Test.sex);
        
        % Preallocate objective values, v, y, and w
        solutionStruct.f = nan(1,length(fullRxnList));
        solutionStruct.y = nan(length(model2Test.mets),length(fullRxnList));
        
        if any(strcmp(model2Test.rxns, strcat('DM_', metabolite{1})))
            solutionStruct.v = nan(length(model2Test.rxns),length(fullRxnList));
            solutionStruct.w = nan(length(model2Test.rxns),length(fullRxnList));
        else
            solutionStruct.v = nan(length(model2Test.rxns)+1,length(fullRxnList));
            solutionStruct.w = nan(length(model2Test.rxns)+1,length(fullRxnList));
        end
        
        % Pre-allocate microbiome biomass coefficients for shadow price
        % analysis
        if any(contains(model2Test.rxns,'Micro_EX_'))
            % Find metabolite coefficients of the community biomass reaction
            communityCoef = full(model2Test.S(:, contains(model2Test.rxns,'communityBiomass')));
            
            % Find the negative coefficients, i.e., the pan taxon biomass
            % metabolites [c] and obtain taxa names.
            solutionStruct.taxonNames = erase(model2Test.mets(communityCoef < 0),'_biomass[c]');
            
            % Find their relative abundance in percentages
            solutionStruct.relAbundances = -(communityCoef(communityCoef < 0));
            
            % Preallocate matrix for species biomass shadow prices
            solutionStruct.shadowPriceBIO = zeros(length(solutionStruct.taxonNames),length(solutionStruct.rxns));
        end
    end
    
    % Ensure all macro metabolites have upper excretion bound of 100000 to
    % prevent infeasibility
    macros = {
        'Diet_EX_energy[d]'
        'Diet_EX_carbohydrate[d]'
        'Diet_EX_protein[d]'
        'Diet_EX_lipid[d]'
        'Diet_EX_sugars[d]'
        'Diet_EX_money[d]'
        };
    
    model2Test = changeRxnBounds(model2Test, macros, 1000000, 'u');
    
    % Set all food exchange reactions to [0,0]
    model2Test = changeRxnBounds(model2Test, model2Test.rxns(contains(model2Test.rxns, 'Food_EX_')), 0, 'b');
    % Accumulate this model's objective values in a local matrix and assign
    % the whole results(i,:,:) slice after the loops. A direct results(i,j,m)
    % write is not a valid sliced variable inside a parfor with nested loops.
    resultSlice = zeros(length(foodItems), length(metabolite));
    % For each diet to be tested
    for j = 1:length(foodItems)
        % Set the food item to its uptake amount, using the database-specific
        % reaction suffix (_usda or _frida from column 2 of the food-items file).
        modelFoodItem = changeRxnBounds(model2Test, strcat('Food_EX_', foodItems(j), '_', foodDb(j)), -foodAmount(j), 'b');
        % Solve the model for each metabolite of interest
        for m = 1:length(metabolite)
            modelPreMadeDietDemand = modelFoodItem;
            % Check if demand reaction is present
            if any(strcmp(modelFoodItem.rxns, strcat('DM_', metabolite{m})))
                if any(strcmp(modelPreMadeDietDemand.rxns, strcat(metabolite{m}, '_synthesis')))
                    modelPreMadeDietDemand = changeRxnBounds(modelPreMadeDietDemand,...
                        modelPreMadeDietDemand.rxns(strcmp(modelPreMadeDietDemand.rxns, strcat(metabolite{m}, '_synthesis'))), 100000, 'u');
                end
            end
            % Change the bound of the demand reaction to 100000
            modelPreMadeDietDemand = changeRxnBounds(modelPreMadeDietDemand, strcat('DM_', metabolite{m}), 100000, 'u');
            % Change the objective to the demand rxn
            modelPreMadeDietDemand = changeObjective(modelPreMadeDietDemand, strcat('DM_', metabolite{m}));
            % Ensure we maximise the FBA
            modelPreMadeDietDemand.osenseStr='max';
            % Solve
            fbaRes = optimizeWBModel(modelPreMadeDietDemand);
            fprintf('m=%d\n',m)
            % Store the FBA result
            if ~isempty(fbaRes.f)
                resultSlice(j,m) = fbaRes.f;
                if saveFullRes
                    solutionStruct.f(1,(j-1)*numel(metabolite)+m)=fbaRes.f;
                end
            else
                resultSlice(j,m) = NaN;
                if saveFullRes
                    solutionStruct.f(1,(j-1)*numel(metabolite)+m)=NaN;
                end
            end
            
            % If required store the full FBA solution
            if saveFullRes
                solutionStruct.solver = fbaRes.solver;
                solutionStruct.osenseStr{1,(j-1)*numel(metabolite)+m} = 'max';
                solutionStruct.stat(1,(j-1)*numel(metabolite)+m)=fbaRes.stat;
                if isfield(fbaRes, 'origStat')
                    solutionStruct.origStat(1,(j-1)*numel(metabolite)+m)={fbaRes.origStat};
                else
                    solutionStruct.origStat(1,(j-1)*numel(metabolite)+m) = "Not available due to solver type";
                end
                
                if fbaRes.stat==1
                    % Save the entire flux vector, all shadow prices, and reduced costs
                    solutionStruct.v(:,(j-1)*numel(metabolite)+m) = fbaRes.v;
                    solutionStruct.y(:,(j-1)*numel(metabolite)+m) = fbaRes.y;
                    solutionStruct.w(:,(j-1)*numel(metabolite)+m) = fbaRes.w;
                    
                    if any(contains(model2Test.rxns,'Micro_EX_'))
                        % Set the names of the biomass metabolites of the
                        % panSpecies models so they can be found
                        pantaxonNamesmassMet = strcat(solutionStruct.taxonNames, '_biomass[c]');
                        % Save species biomass shadow prices
                        solutionStruct.shadowPriceBIO(:,(j-1)*numel(metabolite)+m) = fbaRes.y(matches(model2Test.mets,pantaxonNamesmassMet));
                    end
                end
            end
        end
    end
    % Assign this model's objective values as one sliced write (results(i,:,:)
    % with i the parfor variable and colons elsewhere is a valid slice).
    results(i,:,:) = reshape(resultSlice, 1, length(foodItems), length(metabolite));
    % Save the solution structure. solutionStruct only exists when
    % saveFullRes is true, so guard the save; parsave is used because a
    % bare save() is not permitted inside a parfor loop.
    if saveFullRes
        fname = cell2mat(strcat(resultDir, filesep, 'fbaSol_', wbmModels(i)));
        parsave(fname, solutionStruct);
    end
end

for i = 1:length(metabolite)
    % Extract the results for 1 metabolite from the 3D array
    result2Save = results(:,:,i);
    
    % Set the rows and column names
    rows = wbmModels;
    rows = erase(rows, '.mat');
    
    columns = foodItems;
    columns = strrep(columns, '.txt', '');
    columns = strrep(columns, '.m', '');
    columns = strcat(columns, '_', metabolite{i});
    
    % Convert the array to a table
    result2Save = array2table(result2Save, 'VariableName', columns, 'RowNames', rows);
    % Add sex information
    result2Save.Sex = sex';
    
    % Save the table
    tableFileName = [resultDir, filesep, 'fbaResFood_', metabolite{i}, '.xlsx'];
    writetable(result2Save, tableFileName, 'WriteRowNames', true);
end
end

function parsave(fname, solutionStruct)
% Saves the fields of solutionStruct to fname from within a parfor loop.
% A direct save() call is not transparent and therefore not allowed inside
% parfor, so the save is wrapped in this helper (same approach as
% analyseWBMs.m). The -struct switch saves each field as its own variable.
save(fname, '-struct', 'solutionStruct');
end
