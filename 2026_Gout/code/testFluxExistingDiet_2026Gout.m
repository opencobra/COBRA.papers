function testFluxExistingDiet_2026Gout(diets2Test, metabolite, resultDir, modelPath, saveFullRes, normalise, dietKcals)
% Maximises the production of one or more metabolites on a set of pre-defined
% diets for every whole-body metabolic (WBM) model in a directory, and writes
% the objective values (and optionally the full FBA solutions) to disk.
%
% USAGE:
%
%    testFluxExistingDiet_2026Gout(diets2Test, metabolite, resultDir, modelPath, saveFullRes, normalise, dietKcals)
%
% INPUTS:
%    diets2Test:     Cell array of diet definition file names, each either a
%                    `.m` script or a two-column `.txt` file of diet exchange
%                    reactions and their values.
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
%    normalise:      Logical, if true each diet is scaled to a 2000 kcal
%                    reference using `dietKcals` (default: false).
%    dietKcals:      Numeric vector with one caloric value per diet; required
%                    when `normalise` is true (default: []).
%
% OUTPUTS:
%    None. Files are written to `resultDir`:
%
%                     * fbaSol_<model>.mat - full FBA solution per model (when saveFullRes is true)
%                     * fbaResDiet_<metabolite>.xlsx - objective value per model (row) and diet (column)
%
% NOTE:
%
%    Dietary urate uptake (gut to blood) is not added here; it is added to
%    the food-item WBMs once, by addUrateUptake2WBMs_2026Gout in the driver.
%    Requires the COBRA Toolbox and the WBM/Persephone helper functions
%    (setDietConstraints, addDemandReaction,
%    changeRxnBounds, changeObjective, optimizeWBModel).
%
% .. Author: - Bram Nap

if nargin <  6
    normalise = false;
    dietKcals = [];
end

% When normalising, a caloric value must be supplied for every diet so that
% dietKcals(j) is always defined in the loop below.
if normalise && numel(dietKcals) ~= numel(diets2Test)
    error(['Normalisation requires one caloric value per diet. ' ...
        'There are %d diets but %d caloric values were given.'], ...
        numel(diets2Test), numel(dietKcals));
end

% Find the file names of the wbm models to be tested
wbmModels = what(modelPath).mat;

% Initialise the result variable
results = zeros([length(wbmModels), length(diets2Test), length(metabolite)]);
% Initialise the sex variable
sex = {};

% Define rxn IDs, one per solution column (diet x metabolite)
if isscalar(diets2Test)
    % If only 1 diet is tested the metabolites are unique ids
    fullRxnList = metabolite(:);
else
    % If multiple diets are tested, add diet name to metabolite name for
    % unique rxn ids (<diet>_<metabolite>, same as the result table columns).
    fullRxnList = cell([length(metabolite)*length(diets2Test),1]);
    idxCount = 0;
    for k = 1:length(diets2Test)
        % allocate the names in the list
        startIdx = idxCount*length(metabolite) + 1;
        endIdx = startIdx + length(metabolite) - 1;
        fullRxnList(startIdx:endIdx,1) = strcat(diets2Test(k), '_', metabolite);
        % Increase counter
        idxCount = idxCount + 1;
    end
    % Remove .mat and .txt from the names
    fullRxnList = strrep(fullRxnList, '.m', '');
    fullRxnList = strrep(fullRxnList, '.txt', '');
end

for i = 1:length(wbmModels)
    % Load the model
    model2Test = load(fullfile(modelPath, wbmModels{i}));
    
    % Store the sex of the model
    sex{i} = model2Test.sex;

    % Show which model number is currently being investigated
    disp(i)
    
    % Add the demand reaction
    model2Test = addDemandReaction(model2Test, metabolite);
    model2Test = changeRxnBounds(model2Test, strcat('DM_', metabolite), 0, 'b');
    
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
    % For each diet to be tested
    for j = 1:length(diets2Test)
        % Load the diet based on if it is .m or .txt file
        if contains(diets2Test(j), '.m')
            run(diets2Test{j});
        elseif contains(diets2Test(j), '.txt')
            Diet = readcell(diets2Test{j});
            Diet = Diet(2:end, :);
            Diet(:,2) = cellstr(num2str(cell2mat(Diet(:,2))));
            Diet(:,1) = strcat('Diet_', Diet(:,1));
            Diet(:,1) = strrep(Diet(:,1), '(e)', '[d]');
        end

        % If the diets need to be normalised
        if normalise
            % Calculate the factor required to normalise the met levels
            adjustFactor = dietKcals(j)/2000;
            % Extract met values
            dietValues = Diet(:,2);
            % Convert into int
            dietValues = str2double(string(dietValues));
            % Adjust values, but exempt metabolites below 1 mmol/human/day
            % from scaling (their contribution is negligible and rescaling
            % them only adds noise / precision loss).
            dietValuesAdjuts = dietValues;
            adjustMask = dietValues >= 1 & dietValues <= 100000;
            dietValuesAdjuts(adjustMask) = dietValues(adjustMask)./adjustFactor;
            % convert back into string
            Diet(:,2) = cellstr(num2str(dietValuesAdjuts));
        end

        % Set the dietary constraints on the WBM
        modelPreMadeDiet = setDietConstraints(model2Test, Diet);

        % Solve the model for each metabolite of interest
        for m = 1:length(metabolite)
            % Check if demand reaction is present
            if any(strcmp(modelPreMadeDiet.rxns, strcat('DM_', metabolite{m})))
                modelPreMadeDietDemand = modelPreMadeDiet;
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
                results(i,j,m) = fbaRes.f;
                if saveFullRes
                    solutionStruct.f(1,(j-1)*numel(metabolite)+m)=fbaRes.f;
                end
            else
                results(i,j,m) = NaN;
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
    % Save the solution structure (only exists when saveFullRes is true)
    if saveFullRes
        fname = cell2mat(strcat(resultDir, filesep, 'fbaSol_', wbmModels(i)));
        save(fname,'-struct', 'solutionStruct')
    end
end

for i = 1:length(metabolite)
    % Extract the results for 1 metabolite from the 3D array
    result2Save = results(:,:,i);

    % Set the rows and column names
    rows = wbmModels;
    rows = erase(rows, '.mat');

    columns = diets2Test;
    columns = strrep(columns, '.txt', '');
    columns = strrep(columns, '.m', '');
    columns = strcat(columns, '_', metabolite{i});
    % Convert the array to a table
    result2Save = array2table(result2Save, 'VariableName', columns, 'RowNames', rows);
    % Add sex information
    result2Save.Sex = sex';

    % Save the table
    tableFileName = [resultDir, filesep, 'fbaResDiet_', metabolite{i}, '.xlsx'];
    writetable(result2Save, tableFileName, 'WriteRowNames', true);
end
end
