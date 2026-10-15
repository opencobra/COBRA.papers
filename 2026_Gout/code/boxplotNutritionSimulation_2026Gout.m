function boxplotNutritionSimulation_2026Gout(dataPath, saveDir, fileName, yaxisText, foodMetadataPath, adjust, metaboliteName)
% Draws box plots of predicted metabolite fluxes per diet or food item and
% saves them as `.fig` and `.svg`. Two figures are produced: one grouped by
% sex (female/male) and one without sex separation. When food metadata is
% supplied the x-axis is sorted by food category and coloured bands (Eat /
% Limit / Avoid) are overlaid behind the boxes.
%
% USAGE:
%
%    boxplotNutritionSimulation_2026Gout(dataPath, saveDir, fileName, yaxisText, foodMetadataPath, adjust, metaboliteName)
%
% INPUTS:
%    dataPath:          Char, path to the results table (`.xlsx`); first
%                       column is the model id, the last column is `Sex`, the
%                       columns in between are the diets or food items to plot.
%    saveDir:           Char, directory where the figures are saved.
%    fileName:          Char, base file name; the suffixes `seperatedBySex`
%                       and `noSeperation` are appended.
%    yaxisText:         Char or cell array of char, the y-axis label.
%
% OPTIONAL INPUTS:
%    foodMetadataPath:  Char, path to the food metadata table containing the
%                       `Usda`, `shortname` and `Category` columns. When given
%                       the x-axis is reordered by category and the coloured
%                       bands are drawn (default: '', no adjustment).
%    adjust:            Numeric per-model baseline vector subtracted from the
%                       data before plotting (e.g. the EU average diet flux),
%                       default: '' (no baseline subtraction).
%    metaboliteName:    Char, the metabolite suffix (e.g. `'urate[bc]'`) removed
%                       from the column headers to recover the diet/food label.
%                       Defaults to `'urate[bc]'` for backward compatibility.
%
% OUTPUTS:
%    None. Files written to `saveDir`:
%
%                        * <fileName>seperatedBySex.fig / .svg
%                        * <fileName>noSeperation.fig / .svg
%
% NOTE:
%
%    The x-axis is numeric (positions 1..N) rather than categorical so the
%    boxes and colour bands sit flush against the spines; a categorical ruler
%    keeps an unremovable edge margin.
%
% .. Author: - Bram Nap

% Defaults for the optional inputs.
if nargin < 5 || isempty(foodMetadataPath)
    foodMetadataPath = '';
end
if nargin < 6
    adjust = '';
end
if nargin < 7 || isempty(metaboliteName)
    metaboliteName = 'urate[bc]';
end
adjustXaxis = ~isempty(foodMetadataPath);

% Suffix stripped from the column headers to recover the diet/food label.
metSuffix = ['_' char(metaboliteName)];

% read data
data = readtable(dataPath, "VariableNamingRule","preserve");

% Subtract the per-model baseline (usually EU average diet flux) if supplied.
% Rows are models, columns 2:end-1 are the diets/food items, so the nModel x 1
% baseline broadcasts across the item columns.
if ~isempty(adjust)
    data{:, 2:end-1} = data{:, 2:end-1} - adjust;
end

% Read metadata, only keeping the food items present in the results
if ~isempty(foodMetadataPath)
    foodMetadata = readtable(foodMetadataPath);
    items = strrep(data.Properties.VariableNames(2:end-1), metSuffix, '');
    foodMetadata = foodMetadata(ismember(string(foodMetadata.Usda), items), :);
end

% initialise data storage array
transformedData = cell([2,0]);

% Populate the array per model row so it can be used by boxchart. Every model
% row is included (a boxchart needs the full distribution).
for i = 1:size(data,1)
    % Create new data for row i
    toAdd = [data.Properties.VariableNames(2:end-1)', num2cell(data{i,2:end-1})'];
    toAdd(:,1) = strrep(toAdd(:,1), metSuffix, '');
    % Add sex info
    toAdd(:,3) = data.Sex(i);

    % If data is given, add short names for improved readability of figs
    if adjustXaxis
        [~,~,idx] = intersect(toAdd(:,1), string(foodMetadata.Usda), 'stable');
        toAdd(:,1) = foodMetadata.shortname(idx);
    end
    % Concat
    transformedData = [transformedData; toAdd];
end

% Prepare x tick data
if adjustXaxis
    % Sort the table based on the catogery column indicating how good the
    % food is with 1 being good and 3 being bad
    sortedTable = sortrows(foodMetadata,"Category","ascend");
    % Extract the sorted shortnames
    sorter = sortedTable.shortname;
    % Make the shortnames from the transformedData array a categorical
    % thats ordered based on the sorter we just performed
    xTickInfo = categorical(transformedData(:,1), sorter);

    % x axis limits of the coloured bands (one per category present)
    [bands, colours, bandNames] = getCategoryBands(sortedTable.Category);
else
    % Only transform into caterogical data if no additional info is given
    xTickInfo = categorical(transformedData(:,1));
end

% Match barchartNutritionSimulations_2026Gout.m: use NUMERIC x positions (not a
% categorical ruler, whose unremovable edge margin leaves the colour bands
% short of the spines) and restore the names as tick labels. double(xTickInfo)
% is each observation's position in category order.
xCatLabels = categories(xTickInfo);   % box names, in category order
xPos       = double(xTickInfo);        % numeric x position per observation

% Capitalise sex labels so the legend reads Female / Male.
sexGroups = regexprep(transformedData(:,3), '^.', '${upper($0)}');

f = figure();
hold on
% Make boxplot with sex as seperation
boxchart(xPos, cell2mat(transformedData(:,2)), 'GroupByColor', sexGroups);
% Legend created HERE (before the bands): AutoUpdate appends the band patches
% via their DisplayName, so it reads Female / Male / Eat / Limit / Avoid.
legend()
ylabel(yaxisText)
xticks(1:numel(xCatLabels))
xticklabels(xCatLabels)
xtickangle(90)
xlim([0.5, numel(xCatLabels) + 0.5])   % category edges: bands sit flush

if adjustXaxis
    drawColourBands(bands, colours, bandNames);
end

% Save as svg and fig
savefig(strcat(saveDir, filesep, fileName, 'seperatedBySex.fig'));
saveas(f, strcat(saveDir, filesep, fileName, 'seperatedBySex.svg'));

hold off

f = figure();
hold on
% Make boxplot with no sex separation
bc = boxchart(xPos, cell2mat(transformedData(:,2)));
bc.Annotation.LegendInformation.IconDisplayStyle = 'off';   % keep boxes out of legend
ylabel(yaxisText)
xticks(1:numel(xCatLabels))
xticklabels(xCatLabels)
xtickangle(90)
xlim([0.5, numel(xCatLabels) + 0.5])

if adjustXaxis
    drawColourBands(bands, colours, bandNames);
    legend()                          % show only the band names (boxes excluded above)
end

% Save as svg and fig
savefig(strcat(saveDir, filesep, fileName, 'noSeperation.fig'));
saveas(f, strcat(saveDir, filesep, fileName, 'noSeperation.svg'));

hold off
end

function [bands, colours, bandNames] = getCategoryBands(category)
% x limits, colours and names of the Eat (1) / Limit (2) / Avoid (3) bands
% for a category vector sorted in ascending order. Categories without food
% items are skipped.
allColours = [0.6 1 0.6;   % Green
    1 1 0.6;               % Yellow
    1 0.6 0.6];            % Red
allNames = {'Eat', 'Limit', 'Avoid'};
bands = zeros(0, 2); colours = zeros(0, 3); bandNames = {};
for c = 1:3
    idx = find(category == c);
    if ~isempty(idx)
        bands(end+1, :) = [idx(1) - 0.5, idx(end) + 0.5]; %#ok<AGROW>
        colours(end+1, :) = allColours(c, :); %#ok<AGROW>
        bandNames{end+1} = allNames{c}; %#ok<AGROW>
    end
end
end

function drawColourBands(bands, colours, bandNames)
% Overlay the Eat/Limit/Avoid coloured rectangles behind the current axes. The
% band heights span the current y-limits; call after the data is plotted.
xp = [bands fliplr(bands)];                                   % X-coordinates
yp = ([[1;1]*min(ylim); [1;1]*max(ylim)] * ones(1,size(bands,1))).';   % Y-coordinates
for k = 1:size(bands,1)
    patch(xp(k,:), yp(k,:), colours(k,:), 'FaceAlpha',0.2, ...
        'EdgeColor', colours(k,:), 'DisplayName', bandNames{k});
end
end
