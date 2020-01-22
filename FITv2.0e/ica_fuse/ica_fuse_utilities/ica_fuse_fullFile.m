function files = ica_fuse_fullFile(varargin)
% form full file path using the input directory information and the size of
% the files

inputDir = pwd;
files = [];

if mod(nargin, 2) ~= 0
    error('arguments must be passed in pairs');
end

for ii = 1:2:nargin
    if strcmp(lower(varargin{ii}), 'files')
        files = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'directory')
        inputDir = varargin{ii + 1};
    end
end

temp = repmat(struct('files', []), 1, size(files, 1));
if ~isempty(files)
    for ii = 1:size(files, 1)
        temp(ii).files = fullfile(inputDir, files(ii, :));
    end
    files = str2mat(temp.files); % full file path of the files
else
    files = []; % present working directory
end
        
