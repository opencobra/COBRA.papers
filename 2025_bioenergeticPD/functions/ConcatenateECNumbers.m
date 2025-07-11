function uniqueECString = ConcatenateECNumbers(ecStrings)
    % ecStrings: Cell array where each element contains EC numbers
    
    % Initialize a cell array to store all EC numbers
    allECNumbers = {};
    
    % Regular expression to match EC numbers (complete or incomplete)
    ECpattern = 'EC=\d+\.\d+\.\d+\.\d+|EC=\d+\.\d+\.\d+\.-';
    
    % Loop through each string in the input cell array
    for i = 1:numel(ecStrings)
        % Get the current EC string
        ecString = ecStrings{i};
        
        % Skip if the EC string contains 'EC number not found'
        if contains(ecString, 'EC number not found')
            continue;
        end
        
        % Extract EC numbers from the current string
        matches = regexp(ecString, ECpattern, 'match');
        
        % Remove the 'EC=' prefix from each match
        if ~isempty(matches)
            ecNumbers = regexprep(matches, 'EC=', '');
            allECNumbers = [allECNumbers; ecNumbers(:)]; % Ensure column vector
        end
    end
    
    % Remove duplicates across the entire list
    uniqueECNumbers = unique(allECNumbers);
    
    % Concatenate all unique EC numbers into a single string
    uniqueECString = strjoin(uniqueECNumbers, ', ');
end
