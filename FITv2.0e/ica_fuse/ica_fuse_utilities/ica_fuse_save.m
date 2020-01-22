function ica_fuse_save(varargin)
%% Function to save the MAT files. The function inputs are the same as save command. For versions greater than MATLAB 6.5,
% option is provided to save the MAT files in the appropriate version.
%

%% Load defaults
ica_fuse_defaults;

%% Enforce MAT file version
global ENFORCE_MAT_FILE_VERSION;

if isempty(ENFORCE_MAT_FILE_VERSION)
    ENFORCE_MAT_FILE_VERSION = '-v6';
end

%% MATLAB version
matlab_version = ica_fuse_get_matlab_version;

%% Loop over nargs
for ii = 2:nargin
    if ~strcmp(varargin{ii}(1), '-')
        % get only the variable names
        eval([varargin{ii}, ' = ', 'evalin(''caller'', varargin{ii});']);
    end
end

%% Enforce MAT file version
if (matlab_version > 13)
    save(varargin{1}, varargin{2:end}, ENFORCE_MAT_FILE_VERSION);
else
    save(varargin{1}, varargin{2:end});
end