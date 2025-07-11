function uniqueECList = mergeECNumbers(ecStrings)
    % ecStrings: Cell array of strings, each containing EC numbers
    
    % Initialize a cell array to store all EC numbers
    allECNumbers = {};
    
    % Loop through each string in the input cell array
    for i = 1:numel(ecStrings)
        % Extract EC numbers from the current string
        ecNumbers = extractECNumbers(ecStrings{i});
        
        % Append extracted EC numbers to the list
        if ~isempty(ecNumbers)
            allECNumbers = [allECNumbers; ecNumbers(:)]; % Ensure column vector
        end
    end
    
    % Remove duplicates
    uniqueECNumbers = unique(allECNumbers);
    
    % Create the output cell array
    uniqueECList = cell(length(uniqueECNumbers), 1);
    for i = 1:length(uniqueECNumbers)
        uniqueECList{i} = uniqueECNumbers{i};
    end
end

function ecNumbers = extractECNumbers(data)
    % Extract EC numbers from a single string
    
    % Normalize delimiters by replacing ';', ',' and other potential delimiters
    normalizedData = regexprep(data, '[;,]', ' '); % Replace ';' and ',' with space
    
    % Regular expression to match EC numbers with 'EC:' prefix
    ECpattern = 'EC:\d+\.\d+\.\d+\.\d+';
    
    % Extract matches
    matches = regexp(normalizedData, ECpattern, 'match'); % Find all matches
    
    if isempty(matches)
        ecNumbers = {};
    else
        % Remove the 'EC:' prefix from each match
        ecNumbers = regexprep(matches, 'EC:', '');
        
        % Remove any leading or trailing spaces
        ecNumbers = strtrim(ecNumbers);
    end
end
