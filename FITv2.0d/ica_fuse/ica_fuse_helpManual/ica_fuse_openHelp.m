function ica_fuse_openHelp
% HTML Help for Fusion ICA Toolbox

% get full path for toolbox
toolboxPath = which('ica_fuse_openHelp.m');

% Get the directory where toolbox is running
[pathstr] = fileparts(toolboxPath);

% Source folder where files are located
sourcePath = fullfile(pathstr, 'FITHelpManual.htm');

try
    helpHandle = helpdlg('Opening FIT Help Manual...');

    warning off all;
    s = web(sourcePath, '-browser');
    warning on;
    try
        delete(helpHandle);
    catch
    end
    if s ~= 0
        status = system(['firefox ', sourcePath]);
        if status == 1
            web(sourcePath, '-new');
        end
    end
catch
    disp(lasterr);
end