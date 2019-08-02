function  ica_fuse_check_hdr_file(file_to_check);      
% function to check if there is a header file or not for a image file
%
% Input:
% file_to_check - .img file

[pathstr, fileName, extn] = fileparts(deblank(file_to_check(1, :)));

% Check the header file for analyze images
if strcmpi(extn, '.img')
    
    hdrFile1 = [file_to_check(1:end-3), 'hdr'];
    hdrFile2 = [file_to_check(1:end-3), 'HDR'];
    
    if ~(exist(hdrFile1, 'file') | exist(hdrFile2, 'file'))
        error(['Header file for image ', file_to_check, ' doesn''t exist']);
    end
    
end