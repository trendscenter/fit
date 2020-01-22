function ica_fuse_parallelICA_batch_file(inputFile)
% Batch file for doing parallel ICA fusion


[pathstr, file_name, extn] = fileparts(inputFile);

if ~strcmpi(extn, '.m')
    error(['input file: ', inputFile , ' must be a M-file']);
end

if isempty(pathstr)
    pathstr = fileparts(which(inputFile));
end

inputFile = fullfile(pathstr, [file_name, extn]);

if ~exist(inputFile, 'file')
    error(['input file: ', inputFile , ' doesn''t exist or specify full file path of the file.']);
end

if ~isappdata(0, 'FileRootsData')
    % Use MATLAB Server to get drives on Windows when matlab -nojvm option
    % is used
    ica_fuse_getdrives;
    % Close MATLAB Server
    ica_fuse_closeMatlabServer;
end

% Setup Analysis for parallel ICA
fusionFile = ica_fuse_setup_parallelICA(inputFile);

% Run parallel ICA
ica_fuse_run_parallelICA(fusionFile);

