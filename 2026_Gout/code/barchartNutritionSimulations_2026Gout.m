function barchartNutritionSimulations_2026Gout(dataPath, saveDir, fileName, yaxisText, foodMetadataPath, adjust, metaboliteName)
% Draws a grouped (female/male) bar chart of predicted metabolite fluxes per
% diet or food item and saves it as `.fig` and `.svg`. When food metadata is
% supplied the x-axis is sorted by food category and coloured bands (Eat /
% Limit / Avoid) are overlaid behind the bars.
%
% USAGE:
%
%    barchartNutritionSimulations_2026Gout(dataPath, saveDir, fileName, yaxisText, foodMetadataPath, adjust, metaboliteName)
%
% INPUTS:
%    dataPath:          Char, path to the results table (`.xlsx`); first
%                       column is the model id, the last column is `Sex`, the
%                       columns in between are the diets or food items to plot.
%    saveDir:           Char, directory where the figures are saved.
%    fileName:          Char, base file name; the suffix `seperatedBySex` is
%                       appended.
%    yaxisText:         Char or cell array of char, the y-axis label.
%
% OPTIONAL INPUTS:
%    foodMetadataPath:  Char, path to the food metadata table containing the
%                       `Usda`, `shortname` and `Category` columns. When given
%                       the x-axis is reordered by category and the coloured
%                       bands are drawn (default: '', no adjustment).
%    adjust:            Numeric per-model baseline vector subtracted from the
%                       data before plotting (e.g. the EU average diet flux),
%                       default: ''. When more than two model rows are supplied
%                       and a baseline is given, the rows are collapsed to
%                       female/male means (one bar per sex per item).
%    metaboliteName:    Char, the metabolite suffix (e.g. `'urate[bc]'`) removed
%                       from the column headers to recover the diet/food label.
%                       Defaults to `'urate[bc]'` for backward compatibility.
%
% OUTPUTS:
%    None. Files written to `saveDir`:
%
%                        * <fileName>seperatedBySex.fig
%                        * <fileName>seperatedBySex.svg
%
% NOTE:
%
%    The x-axis is numeric (positions 1..N) rather than categorical so the bars
%    and colour bands sit flush against the spines; a categorical ruler keeps an
%    unremovable edge margin.
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

% Adjust data based on input, which is usually baseline eu average values.
% Rows are models, columns 2:end-1 are the diets/food items, so the nModel x 1
% baseline broadcasts across the item columns.
if ~isempty(adjust)
    data{:, 2:end-1} = data{:, 2:end-1} - adjust;
    if size(data,1) > 2
        % Collapse many (microbiome) models into female/male means so the bar
        % chart shows one bar per sex per item. Keep the first two rows to
        % preserve the column schema, then overwrite them with the means.
        avgFemale = mean(data{strcmp(data.Sex,'female'), 2:end-1}, 1);
        avgMale   = mean(data{strcmp(data.Sex,'male'),   2:end-1}, 1);
        data = data(1:2, :);
        data{1, 2:end-1} = avgFemale;
        data{2, 2:end-1} = avgMale;
        data.Sex = {'Female'; 'Male'};
    end
end
% Read metadata
if ~isempty(foodMetadataPath)
    foodMetadata = readtable(foodMetadataPath);
end

% Barchart requires data to be a column for each x-axis and row for each
% additional bar per x-axis tick

% Set the x ticks and format: strip file extensions and the metabolite suffix.
xTicks = strrep(data.Properties.VariableNames(2:end-1), '.m', '');
xTicks = strrep(xTicks, '.txt', '');
xTicks = strrep(xTicks, metSuffix, '');

if adjustXaxis
    % Only keep the metadata of the food items present in the results
    foodMetadata = foodMetadata(ismember(string(foodMetadata.Usda), xTicks), :);
    % Sort the table based on the catogery column indicating how good the
    % food is with 1 being good and 3 being bad
    sortedTable = sortrows(foodMetadata,"Category","ascend");
    % Extract the sorted shortnames
    sorter = sortedTable.shortname;

    % Find the shortNames
    [~, ~, idx] = intersect(xTicks, string(sortedTable.Usda), 'stable');

    % Replace old headers with the shortnames
    data.Properties.VariableNames(2:end-1) = sortedTable.shortname(idx);
    dataSorted = data(:, sorter);

    % x axis limits of the coloured bands (one per category present)
    [bands, colours, bandNames] = getCategoryBands(sortedTable.Category);
    xLabels = dataSorted.Properties.VariableNames;   % names, already in sort order
else
    xLabels = xTicks;
    dataSorted = data(:, 2:end-1);
end

f = figure();
hold on
% Plot on a NUMERIC x-axis (positions 1..N) instead of a categorical one. A
% categorical ruler always reserves an edge margin that cannot be removed
% (XLimitMethod is a no-op on it), leaving white space between the end bars and
% the spines. A numeric axis lets us set the limits so the bars and colour
% bands sit flush. The names are restored as tick labels afterwards.
xPos = 1:numel(xLabels);
bar(xPos, dataSorted{:, :}')        % N groups (one per item) x 2 series (sex)
% Create the legend HERE (right after bar), not at the end: creating it after
% xtickangle/xlim changes the saved axes state and makes the figure overflow
% when reparented into the tiledlayout. Sex labels get a capitalised first
% letter; AutoUpdate (on by default) appends the band patches below via their
% DisplayName, so they show as Safe/Medium/Harmful instead of data1/2/3.
sexNames = regexprep(cellstr(string(data.Sex)), '^.', '${upper($0)}');
legend(sexNames)
ylabel(yaxisText)
xticks(xPos)
xticklabels(xLabels)
xtickangle(90)
xlim([0.5, numel(xPos)+0.5])        % category edges: bars sit flush, no extra margin
% ylim([400,525]);
if adjustXaxis
    drawColourBands(bands, colours, bandNames);
end

% Save as svg and fig
savefig(strcat(saveDir, filesep, fileName, 'seperatedBySex.fig'));
saveas(f, strcat(saveDir, filesep, fileName, 'seperatedBySex.svg'));

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
