function [data, V] = ica_fuse_loadData(P, file_numbers)
% loads the data for .img, .nii, ascii files
%
% Input:
% 1. P - file names
% 2. file_numbers - File Numbers
% Output
% data - can be a spatial map or a time course
%
% Image data will be of dimensions xdim, ydim, zdim, N (number of files)
% Ascii data will be of dimensions data length, 2, N


firstFile = deblank(P(1, :));

[firstFile, number] = ica_fuse_parseExtn(firstFile);

% check the extension of the first file
[pathstr, fileN, extn] = fileparts(firstFile);

if strcmpi(extn, '.img')
    ica_fuse_check_hdr_file(firstFile);
end

if ~exist('file_numbers', 'var')
    file_numbers = [];
end


if strcmpi(extn, '.img') | strcmpi(extn, '.nii')
    
    % get the volume information
    V = ica_fuse_spm_vol(P);
    if ~isempty(file_numbers)
        V = V(file_numbers);
    end
    
    % read data
    data = ica_fuse_read_vols(V);
    
else
    
    if ~isempty(file_numbers)
        P = P(file_numbers, :);
    end
    
    % treat file as an ascii file
    for ii = 1:size(P, 1)
        
        currentFile = deblank(P(ii, :));
        % load ascii file
        temp = load('-ascii', currentFile);
        
        % Make data a column vector
        if prod(size(temp)) == length(temp)
            if size(temp, 1) == 1
                temp = temp';
            end
        end
        % end for making data a column vector
        
        % First column will be time or scans followed by data
        if size(temp, 2) > 1
            dat = temp(:, 1:2);
        else
            % add scans as first column
            dat = [(1:1:length(temp))', temp];
        end
        
        clear temp;
        
        size_dat = size(dat);
        
        % check the dimensions of each file with that of first file
        if ii == 1
            check_size = size_dat;
            data = zeros(size(dat, 1), 2, size(P, 1));
        else
            if length(find((check_size == size_dat) > 0)) ~= length(size_dat)
                error(['Dimensions of file ', currentFile, ' doesn''t match that of first file ', deblank(P(1, :))]);
            end
        end
        
        data(:, :, ii) = dat;
        clear dat;
    end
    % end loop over files
    
    V(1).dim = size(data, 1);
    
end
% end for loading data

% replace NaN with zero
data(isnan(data)) = 0;