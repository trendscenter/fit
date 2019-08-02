function ica_fuse_printString(fid, titleToPrint, stringToPrint)
% function to print string to screen or file
% 
% Input:
% 1. fid - file id
% 2. titleToPrint - Title string to print
% 3. stringToPrint - String to print
% 

if ~isempty(titleToPrint)
    fprintf(fid, '%s\n', titleToPrint);
end

for ii = 1:size(stringToPrint, 1)
    fprintf(fid, '%s\n', deblank(stringToPrint(ii, :)));
end
