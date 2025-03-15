function ECnumbers = getECnumbersFromUniProt(entryIDs)
    % Check if the input is a cell array of strings or a string array
    if ~iscell(entryIDs) && ~isstring(entryIDs)
        error('Entry IDs must be a cell array of strings or a string array.');
    end

    % Initialize the output as a cell array
    numEntries = length(entryIDs);
    ECnumbers = cell(numEntries, 1);

    % UniProt base URL for querying
    baseURL = 'https://www.uniprot.org/uniprot/';

    % Loop over each Entry ID
    for i = 1:numEntries
        entryID = char(entryIDs(i));  % Get the current Entry ID
        queryURL = [baseURL, entryID, '.txt'];  % Construct the full URL

        % Use Java to fetch the data with proper encoding
        try
            % Create a URL object
            url = java.net.URL(queryURL);

            % Open a connection to the URL
            connection = url.openConnection();
            
            % Set up the input stream and reader with UTF-8 encoding
            inputStream = connection.getInputStream();
            inputStreamReader = java.io.InputStreamReader(inputStream, 'UTF-8');
            bufferedReader = java.io.BufferedReader(inputStreamReader);

            % Read the content line by line
            data = '';
            line = bufferedReader.readLine();
            while ~isempty(line)
                data = [data, char(line), newline];
                line = bufferedReader.readLine();  % Read next line
            end

            % Close the streams
            bufferedReader.close();
            inputStreamReader.close();
            inputStream.close();
        catch
            warning(['Failed to retrieve data for Entry ID: ', entryID, '. Check the Entry ID and your internet connection.']);
            ECnumbers{i} = 'Error retrieving data';
            continue;
        end

        % Search for all EC numbers in the retrieved data
        ECnumbers{i} = extractECnumbers(data);
    end
end

function ECList = extractECnumbers(data)
    % Extract EC numbers from the UniProt data
    ECpattern = 'EC=';
    matches = regexp(data, [ECpattern, '[^\s;]*'], 'match'); % Find all matches

    % Remove duplicates
    if isempty(matches)
        ECList = 'EC number not found';
    else
        uniqueECs = unique(matches); % Remove duplicates
        ECList = strjoin(uniqueECs, ', '); % Join unique EC numbers into a single string
    end
end
