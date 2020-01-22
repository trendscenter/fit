function inputData = ica_fuse_read_variables(file_name, keywd, varCriteria)
% reads variable or variables from file with the specified keywd

% get the file name and extensions
[pathstr, fName, extns] = fileparts(file_name);
    
% initialise input data structure
inputData = struct;

% convert the variable to cell format
if ischar(keywd)
    keywd = cellstr(keywd);
end

% convert the variable to cell format
if ischar(varCriteria)
    varCriteria = cellstr(varCriteria);
end


% check the input file extension
if ~strcmpi(extns, '.m')
    error(['Input file should be an M-file']);
else
    % old directory
    oldDir = pwd;
    if ~isempty(pathstr)
        % change to new directory
        cd(pathstr);    
    end
    inputData = struct;
    % evaluate the input text
    eval(fName);
    % change to old directory
    cd(oldDir);
    % loop over required keywords
    for ii = 1:length(keywd)
        
        currentVar = {};
        try
            currentVar = eval(keywd{ii});
        catch
            error(['Variable ', keywd{ii}, ' doesn''t exist or is not of the required data type in file ', file_name]);
        end  
        
        
        checkVar = currentVar;
        
        % check if is a valid integer or not
        if isnumeric(checkVar)
            checkVar = num2str(checkVar);                        
        end
        
        % get the status and message                        
        [status, message] = ica_fuse_errorCheck(checkVar, varCriteria{ii});
        if status == 0
            error(message);
        end
        % end for getting the correct data type
        inputData = setfield(inputData, keywd{ii}, currentVar);
        clear currentVar; clear checkVar;
    end
    % end for loop over required keywords
end