function saveCorrTableRestore_2026Gout(maxCorr, minCorr, corrSaveDir)
% Save the best positive (Max) and negative (Min) cluster correlations per
% reaction, as obtained from correlateClustersRestore_2026Gout, to a single
% clusterCorrelations.xlsx file. Every reaction/direction is written as its
% own block (header row + one row per cluster size), separated by an empty
% row. Directions without any correctly signed cluster are not written.
%
% How to use:
%   saveCorrTableRestore_2026Gout(maxCorr, minCorr, corrSaveDir)
% Inputs:
%   maxCorr:     Cell array, bestCorrPerSizeList output of
%                correlateClustersRestore_2026Gout
%   minCorr:     Cell array, minbestCorrPerSizeList output of
%                correlateClustersRestore_2026Gout
%   corrSaveDir: Text, directory where clusterCorrelations.xlsx is saved
% Authors:
%   Bram Nap: Feb 2026

% Collect the blocks to save
blocks = {};
for i = 1:size(maxCorr,2)
    for corrTable = {maxCorr{i}, minCorr{i}}
        if istable(corrTable{1}) && height(corrTable{1}) > 0
            blocks{end+1} = [corrTable{1}.Properties.VariableNames; table2cell(corrTable{1})]; %#ok<AGROW>
        end
    end
end

% Pad all blocks to the same number of columns and stack them with an
% empty row in between
nCols = max([cellfun(@(x) size(x,2), blocks), 1]);
cell2Save = {};
for i = 1:numel(blocks)
    block = blocks{i};
    block(:, end+1:nCols) = {' '};
    if ~isempty(cell2Save)
        cell2Save = [cell2Save; repmat({' '}, 1, nCols)]; %#ok<AGROW>
    end
    cell2Save = [cell2Save; block]; %#ok<AGROW>
end

% Save the cell array as xlsx file. Delete an old file first, as writecell
% does not clear cells outside the newly written range.
savePath = strcat(corrSaveDir, filesep, 'clusterCorrelations.xlsx');
if isfile(savePath)
    delete(savePath);
end
writecell(cell2Save, savePath);
end
