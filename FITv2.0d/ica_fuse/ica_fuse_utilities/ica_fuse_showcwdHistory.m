function [previousDir] = ica_fuse_showcwdHistory

% get the recent previous directories. previousDir contains pwd, GIFT path, and the
% other recent directories

% Pull the previous working directories in Matlab by reading the file
% cwdhistory.m

try
    
    % get the old directory
    oldDir = deblank(pwd);
    
    % Fusion path
    % get the path for the Fusion Toolbox
    whichFusion = which('fusion.m');
    fusionPath = fileparts(whichFusion);
    fusionPath = deblank(fusionPath);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get the preferences directory
    preferenecesDir = prefdir;
    
    preferenecesDir = deblank(preferenecesDir);
    
    % Contents of the folder in which preferences are present
    [fileContents] = ica_fuse_listFiles_inDir(preferenecesDir, '*.m');
    
    checkCwdHistory = strmatch('cwdhistory.m', lower(fileContents), 'exact');
    
    % If the cwdhistory file is not present
    % return current working directory and the fusion Path
    if isempty(checkCwdHistory)
        previousDir{1} = oldDir;
        previousDir{2} = fusionPath;
    else
        % open current working directory history m file and read the contents
        fid = fopen(fullfile(preferenecesDir, 'cwdhistory.m'), 'r');
        temp = {};
        if fid ~= -1
            % Read the contents of the cwdhistory and
            [temp, nLines, nonSpace] = ica_fuse_readContents_file(fid);
            fclose(fid);
        end
        
        if isempty(temp)
            error('No directories exist');
        end
        
        
        % Pull only five recent working directories
        if nonSpace > 5
            for ii = 1:5
                previousDir{ii} =  temp{ii};
            end
            % Check if the fusion path ia already present
            checkfusionPath = strmatch(lower(fusionPath), str2mat(lower(previousDir)), 'exact');
            % If fusion path is not present
            if isempty(checkfusionPath)
                previousDir{6} = fusionPath;
            end
        else
            previousDir = temp;
            if isempty(checkfusionPath)
                previousDir{length(temp) + 1} = fusionPath;
            end
        end
    end
    
    % preserve the old directory
    cd(oldDir);
    
    % Check the present working directory in the previous directories
    checkPwd = strmatch(lower(oldDir), str2mat(lower(previousDir)), 'exact');
    
    % Add present working directory
    if isempty(checkPwd)
        previousDir{length(previousDir) + 1} = oldDir;
    end
    
    count = length(previousDir);
    
catch
    
    previousDir{1} = oldDir;
    previousDir{2} = fusionPath;
    count = length(previousDir);
    
    % preserve the old directory
    cd(oldDir);
end
