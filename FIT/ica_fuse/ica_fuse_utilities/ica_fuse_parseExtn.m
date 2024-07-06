function [file, number] = ica_fuse_parseExtn(file)
% Parse extension

[pathstr, fName, extn] = fileparts(file);

commaPos = find(extn == ',');

number = 1;

if ~isempty(commaPos)
    number   = str2num(extn((commaPos(1)+1):end));
    extn = extn(1:(commaPos-1));
    file   = fullfile(pathstr, [fName, extn]);
end
