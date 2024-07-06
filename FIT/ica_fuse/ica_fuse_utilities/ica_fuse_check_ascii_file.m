function ica_fuse_check_ascii_file(currentFile)
% Check if the file input is an ascii or not

try
    load(currentFile, '-ascii');
catch
    error('Error:Ascii', ['Error in reading the file: %s. \nCheck if it is an ascii file.'], ...
        currentFile);
end