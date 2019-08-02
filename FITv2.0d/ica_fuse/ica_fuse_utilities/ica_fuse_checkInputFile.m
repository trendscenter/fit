function [inputText] = ica_fuse_checkInputFile(inputFile, inputText)
% Read from the input file the variables that match the tag in the input
% text and set the answer string for edit or value for the popup controls.
%
% Input:
% 1. inputFile - M file containing the necessary variables to run the batch
% script
% 2. inputText - inputText is the structure containing the information
% regarding user interface controls and the values
%
% Output:
% inputText - Updated inputText in the field answerString for edit or value
% for popup

[pathstr, inputF, extn] = fileparts(inputFile);

if ~strcmpi(extn, '.m')
    error(['Input file: ', inputFile, ' must be a M file']);
end

oldDir = pwd;

try

    if ~isempty(pathstr)
        cd(pathstr);
    end

    % Evaluate the input file
    eval(inputF);

    controlsTobeRead = find([inputText.read_from_file] == 1);

    % Loop over the required controls
    for ii = 1:length(controlsTobeRead)
        % current tag
        currentTag = inputText(controlsTobeRead(ii)).tag;
        dataType  = inputText(controlsTobeRead(ii)).answerType;
        uiType = inputText(controlsTobeRead(ii)).uiType;
        answerString = str2mat(inputText(controlsTobeRead(ii)).answerString);
        try
            currentData = eval(currentTag);
        catch
            error(['Error:', currentTag], ['Check input file: %s\nto make sure that variable %s exists or is of the required data type'], ...
                inputFile, currentTag);
        end

        % For edit controls
        if strcmpi(uiType, 'edit')
            % convert to string
            if isnumeric(currentData)
                currentData = num2str(currentData);
            end
            inputText(controlsTobeRead(ii)).answerString = currentData;

        elseif strcmpi(uiType, 'popup') | strcmpi(uiType, 'listbox')

            % Set the value for popup control
            if strcmp(currentTag, 'maskFile')
                if ~isempty(currentData)
                    if ~strcmp(currentData, '[]')
                        inputText(controlsTobeRead(ii)).value = 2;
                    else
                        inputText(controlsTobeRead(ii)).value = 1;
                    end
                else
                    inputText(controlsTobeRead(ii)).value = 1;
                end

            else

                % Check the data type for the pop up control
                if isnumeric(currentData)
                    if size(answerString, 1) < currentData
                        error(['Error:', currentTag], ['Value for the variable %s in the input file:\n %s exceeds the maximum'], ...
                            currentTag, inputFile);
                    elseif currentData == 0
                        inputText(controlsTobeRead(ii)).value = 1;
                    else
                        inputText(controlsTobeRead(ii)).value = round(currentData);
                    end
                else
                    % for popup set the value
                    matchIndex = strmatch(lower(currentData), lower(answerString), 'exact');
                    if ~isempty(matchIndex)
                        inputText(controlsTobeRead(ii)).value = matchIndex;
                    else
                        if strcmp(currentTag, 'z_scores')
                            matchIndex = strmatch(lower(currentData), {'no', 'yes'}, 'exact');
                            inputText(controlsTobeRead(ii)).value = matchIndex;
                        end
                    end

                end
                % end for checking the data type for the pop up control
            end
            % end for setting the value for popup control
        end
        % end for checking controls

    end
    % End loop over the required controls

    cd(oldDir);

catch

    cd(oldDir);
    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;
end