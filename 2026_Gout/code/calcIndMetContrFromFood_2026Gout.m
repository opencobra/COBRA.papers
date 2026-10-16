function [result] = calcIndMetContrFromFood_2026Gout(foodItemWbm, dietMetList, metabolite)

% Load the model
model = load(foodItemWbm);

if isfield(model, 'model')
    model = model.model;
end

% Set food exchange to 0
model.lb(contains(model.rxns,'Food_EX_'))=0;
model.ub(contains(model.rxns,'Food_EX_'))=0;

% Set the EU Average diet
model = setDietConstraints(model, 'EUAverageDiet');

% Ensure macros exchanges are open to prevent infeasibility
macros = {
    'Diet_EX_energy[d]'
    'Diet_EX_carbohydrate[d]'
    'Diet_EX_protein[d]'
    'Diet_EX_lipid[d]'
    'Diet_EX_sugars[d]'
    'Diet_EX_money[d]'
    };

model = changeRxnBounds(model, macros, 1000000, 'u');

% Dietary urate uptake (gut to blood) is already added to the food-item
% WBMs by addUrateUptake2WBMs_2026Gout in the driver

% Copy the original model
modelEU = model;

% Set up the model to add demand reaction for metabolite of interest
modelEU = addDemandReaction(modelEU, metabolite);
modelEU = changeRxnBounds(modelEU, strcat('DM_' , metabolite), 100000, 'u');
modelEU = changeObjective(modelEU, strcat('DM_' , metabolite));
modelEU.osenseStr = 'max';


cnt=1;
% For every food item to test
for j = 1:length(dietMetList)
    % get the breakdown reaction
    rxnName = strcat('Diet_EX_', dietMetList{j}, '[d]');

    % Remove the macro metabolites from the met list as they not play a role in
    % metabolites and therfore are wasted simulations.
    % macros = strrep(macros, 'Diet_EX_', '');

    % [~, idx, ~] = intersect(mets, macros, 'stable');
    % mets(idx) = [];
    % Initialise a counter

    % Copy model for testing indivudual metabolites
    modelTest = modelEU;

    % Obtain the index for the dietary exchange reaction of
    % metabolite i
    rxnId = find(strcmp(modelTest.rxns,rxnName));
    % add 100 mmol/day of the metabolite on top of the EU diet
    modelTest.lb(rxnId) = -100 + modelTest.lb(rxnId);
    modelTest.ub(rxnId) = 0;
    % Solve model and store results
    fba = optimizeWBModel(modelTest);
    result{cnt,1} = dietMetList{j};
    if fba.f > 0
        result{cnt,2} = fba.f;
    else
        result{cnt,2} = 0;
    end
    % Increase counter
    cnt = cnt +1;

    % Store the cell array in the result structure. X has to be placed in
    % front of the food ID as MATLAB does not allow only numbers as
    % fieldname
    % resultStruct.Individual.(strcat('x',string(dietMetList(j)))) = result;

    % % Store results of all the foods simulations in one table
    % if j == 1
    %     combinedRes = cell2table(result(:, [1 3]));
    %     combinedRes.Properties.VariableNames = {'metID', char(dietMetList(j))};
    %     combinedRes{:,2} = num2cell(str2num(char(combinedRes{:,2})));
    % else
    %     combinedRes2Add = cell2table(result(:, [1 3]));
    %     combinedRes2Add.Properties.VariableNames = {'metID', '1234'};
    %     combinedRes2Add{:,2} = num2cell(str2num(char(combinedRes2Add{:,2})));
    %
    %     combinedRes = outerjoin(combinedRes,combinedRes2Add, 'Keys', 'metID','MergeKeys',true);
    % end
end
% Save the end result in the result structure array
% resultStruct.combinedRes = combinedRes;
end