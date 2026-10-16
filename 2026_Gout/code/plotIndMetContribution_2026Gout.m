function plotIndMetContribution_2026Gout(gfFilePath, mwbmFilePath, saveDir, fileName, highlightMets)
% Plots the increase in urate production after supplementing 100 mmol of
% each individual dietary metabolite (results of calcIndMetContrFromFood_2026Gout,
% driver step 6) as heatmaps and saves them as `.fig`, `.png` and `.svg`.
%
% Panel A: increase per metabolite for the germ-free female and male WBM and
% the average female and male mWBM. Panel B: increase per mWBM for a few
% selected metabolites.
%
% USAGE:
%
%    plotIndMetContribution_2026Gout(gfFilePath, mwbmFilePath, saveDir, fileName, highlightMets)
%
% INPUTS:
%    gfFilePath:      Char, dietMetaboliteContributions.xlsx of the germ-free
%                     WBMs (first column MetID, then one column per model).
%    mwbmFilePath:    Char, dietMetaboliteContributions.xlsx of the mWBMs.
%    saveDir:         Char, directory where the figure is saved.
%    fileName:        Char, base file name of the figure.
%
% OPTIONAL INPUTS:
%    highlightMets:   Cell array of VMH metabolite IDs shown per mWBM in
%                     panel B (default: {}, panel B is not drawn).
%
% NOTE:
%
%    The files hold the maximal urate flux with each metabolite added. The
%    increase is the flux minus the baseline of that model, i.e. the flux
%    without any metabolite added. The baseline is the lowest value of the
%    model's column: adding a metabolite cannot lower the maximal urate flux,
%    and metabolites without effect give exactly the baseline (checked on the
%    full results: lowest value = most frequent value for every model). Increases below 1e-5 are set to 0 and metabolites without an
%    increase in any panel A column are not shown. Model columns are
%    assigned to female or male by '_female' / '_male' (or 'female' / 'male')
%    in the column name. Metabolite names are taken from the VMH database
%    (loadVMHDatabase). Requires MATLAB R2023a or later (sky colormap).
%
% .. Author: - Bram Nap

if nargin < 5
    highlightMets = {};
end
fSize = 9;

% Increases per model
[gfIds, gfInc, gfNames] = readIncrease(gfFilePath);
[mIds, mInc, mNames] = readIncrease(mwbmFilePath);

% Panel A: germ-free female/male and average mWBM female/male
isFem = @(n) contains(n, 'female', 'IgnoreCase', true);
isMale = @(n) contains(n, 'male', 'IgnoreCase', true) & ~isFem(n);
allIds = union(gfIds, mIds, 'stable');
panelA = zeros(numel(allIds), 4);
[~, ia, ib] = intersect(allIds, gfIds, 'stable');
panelA(ia, 1) = mean(gfInc(ib, isFem(gfNames)), 2);
panelA(ia, 2) = mean(gfInc(ib, isMale(gfNames)), 2);
[~, ia, ib] = intersect(allIds, mIds, 'stable');
panelA(ia, 3) = mean(mInc(ib, isFem(mNames)), 2);
panelA(ia, 4) = mean(mInc(ib, isMale(mNames)), 2);
panelA(abs(panelA) < 1e-5 | isnan(panelA)) = 0;

keep = any(panelA ~= 0, 2);
ids = allIds(keep);
panelA = panelA(keep, :);
names = metNames(ids);
[names, order] = sort(names);
panelA = panelA(order, :);

% Panel B: selected metabolites per mWBM
highlightMets = highlightMets(ismember(highlightMets, mIds));
drawB = ~isempty(highlightMets);

fig = figure('Units', 'centimeters', 'Position', [8.11, 1.52, 24, 23.7]);
if drawB
    % 4x2 layout: panel A spans 3 rows x 2 columns, panel B 1 row x 2 columns
    t = tiledlayout(4, 2);
    t.TileSpacing = 'tight';
    t.Padding = 'tight';
    % Leave a strip at the top for the panel label of panel A
    t.Units = 'normalized';
    t.OuterPosition = [0 0 1 0.95];
    nexttile(1, [3 2]);
end
hm1 = heatmap({'GF Female', 'GF Male', 'Avg microbiome Female', 'Avg microbiome Male'}, ...
    names, panelA, 'ColorLimits', [0, 125], 'Colormap', sky(15));
hm1.Title = 'Increase in urate production after supplementation with dietary metabolites';
hm1.FontSize = fSize;

if drawB
    [~, idx] = ismember(highlightMets, mIds);
    modelLabels = erase(mNames, {'mWBM_', 'Day_0_P1', 'foodItem'});
    modelLabels = strrep(modelLabels, '_female', ' F');
    modelLabels = strrep(modelLabels, '_male', ' M');
    modelLabels = strtrim(strrep(modelLabels, '_', ' '));
    nexttile(7, [1 2]);
    hm2 = heatmap(modelLabels, metNames(highlightMets), mInc(idx, :), ...
        'ColorLimits', [0, 25], 'Colormap', sky(10));
    numWords = {'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine'};
    nMets = numel(highlightMets);
    if nMets <= numel(numWords); nText = numWords{nMets}; else; nText = num2str(nMets); end
    if nMets == 1; metText = 'metabolite'; else; metText = 'metabolites'; end
    hm2.Title = sprintf('Individual increase in urate production after supplementation with %s dietary %s', nText, metText);
    hm2.FontSize = fSize;
    addHeatmapLabel(hm1, 'A', fSize);
    addHeatmapLabel(hm2, 'B', fSize);
else
    warning('plotIndMetContribution_2026Gout:noHighlight', ...
        'None of the highlighted metabolites are in the mWBM results; panel B is not drawn.');
end

if ~isfolder(saveDir)
    mkdir(saveDir);
end
savefig(fig, fullfile(saveDir, [fileName, '.fig']));
saveas(fig, fullfile(saveDir, [fileName, '.png']));
saveas(fig, fullfile(saveDir, [fileName, '.svg']));
end

function [ids, increase, modelNames] = readIncrease(filePath)
% Read a dietMetaboliteContributions file and return the increase in urate
% flux of every metabolite (rows) for every model (columns).
data = readcell(filePath);
header = string(data(1, :));
data = data(2:end, :);
ids = cellstr(string(data(:, 1)));
% Model columns: names containing female or male (other columns are ignored)
modelCols = find(contains(header, 'male', 'IgnoreCase', true));
modelNames = cellstr(header(modelCols));
values = cellfun(@toNumber, data(:, modelCols));
% Baseline per model: the lowest value. Adding a dietary metabolite only
% gives the model more to take up, so it cannot lower the maximal urate
% flux; metabolites without effect are at the baseline.
baseline = min(values, [], 1);
increase = values - baseline;
end

function x = toNumber(v)
% Convert a readcell value to a double (missing/text -> NaN)
if isnumeric(v) && isscalar(v)
    x = double(v);
else
    x = str2double(string(v));
end
end

function names = metNames(ids)
% Readable metabolite names from the VMH database; IDs not found are kept
mets = loadVMHDatabase().metabolites;
names = ids(:);
[found, idx] = ismember(ids, mets(:, 1));
names(found) = mets(idx(found), 2);
names = strrep(names, 'tetrahydrofolate', '-THF');
names = strrep(names, '(2-)', '');
names = strrep(names, '(1-)', '');
names = strrep(names, '(4+)', '');
names = regexprep(strtrim(names), '^.', '${upper($0)}');
end

function addHeatmapLabel(hm, lbl, fSize)
% Bold panel label just above the top-left corner of a heatmap. Heatmaps
% are charts (not axes), so a figure annotation at the rendered pixel
% position is used.
fig = ancestor(hm, 'figure');
drawnow;
hmPx = getpixelposition(hm, true);
figPx = getpixelposition(fig);
x = hmPx(1) / figPx(3);
y = (hmPx(2) + hmPx(4)) / figPx(4);
w = 0.05; h = 0.05;
annotation(fig, 'textbox', [min(max(x, 0), 1 - w), min(max(y, 0), 1 - h), w, h], ...
    'String', lbl, 'FontWeight', 'bold', 'FontSize', fSize + 3, 'EdgeColor', 'none', ...
    'Margin', 0, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
