function formattedStr = ica_fuse_formatStr(names, characterToPut)
%% Format String
%
% Inputs:
% 1. names - Cell array or character array of names
% 2. characterToPut - Character you want to put
%
% Outputs:
% formattedStr - Formatted string
%

if (~exist('characterToPut', 'var'))
    characterToPut = ', ';
end

names = cellstr(names);
formatStr = repmat(['%s', characterToPut], 1, length(names));
formattedStr = sprintf(formatStr, names{:});

% Find the last position
pos = ica_fuse_findstr(formattedStr, characterToPut);
formattedStr = formattedStr(1:pos(end)-1);