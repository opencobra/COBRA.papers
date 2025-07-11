function [IDList] = GeneNameToEntrezID(GeneNameList)
% Convert a list of Gene Names to Entrez IDs
% Input: GeneNameList - cell array of gene names
% Output: IDList - cell array of Entrez IDs

% Ensure input is a cell array
if iscell(GeneNameList)
    geneNames = GeneNameList;
else
    geneNames = cellstr(GeneNameList);
end

% Initialize output cell array
IDList = cell(length(geneNames), 1);

% Web options with infinite timeout
options = weboptions('Timeout', Inf);

% Loop through each gene name
for i = 1:length(geneNames)
    geneName = geneNames{i};
    
    if ~ismissing(geneName)
        entrezID = getEntrezID(geneName, options);
        IDList{i} = entrezID;
    end
end
end

function entrezID = getEntrezID(geneName, options)
% Initialize output
entrezID = '';

% Base URLs for Ensembl and NCBI
ensemblURL = ['https://rest.ensembl.org/xrefs/symbol/homo_sapiens/', geneName, '?content-type=application/json'];
ncbiURL = ['https://www.ncbi.nlm.nih.gov/gene/?term=', geneName];

% Try to fetch from Ensembl
try
    response = webread(ensemblURL, options);
    if isstruct(response)
        for k = 1:length(response)
            if strcmp(response(k).dbname, 'EntrezGene')
                entrezID = response(k).primary_id;
                return;
            end
        end
    end
catch
    % If Ensembl fails, proceed to NCBI
    entrezID =[];
end

% Try to fetch from NCBI if Ensembl didn't work
if isempty(entrezID)
    try
        response = webread(ncbiURL, options);
        idStart = strfind(response, 'Gene ID:');
        if ~isempty(idStart)
            % Adjusted to extract only the numeric part after "Gene ID: "
            numericPart = regexp(response(idStart:end), '\d+', 'match');
            if ~isempty(numericPart)
                entrezID = numericPart{1};  % Extract the first matched numeric value
            end
        end
    catch
        % Handle any failure cases (e.g., webread failure)
        entrezID = [];
    end
end

% If both fail, return an empty string or some default value
end