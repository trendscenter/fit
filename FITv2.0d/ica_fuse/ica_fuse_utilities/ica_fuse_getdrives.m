function [d, errorMessage] = ica_fuse_getdrives
% Uses java to return the drives in Windows
% Don't use exist(d, 'dir') as it may take long time

errorMessage = '';

if isappdata(0, 'FileRootsData')
    d = getappdata(0, 'FileRootsData');
    return;
end

% Open matlab server for -nojvm
ica_fuse_openMatlabServer;
global COMSERVER_HANDLE;

if ispc
    try
        % Use java to get drives
        commandToExecute = ['import java.io.*; f = File(''''); a = f.listRoots; d = cell(length(a), 1);', ...
            'for nD = 1:length(d); d{nD} = char(a(nD).toString); end;'];
        if usejava('jvm')
            eval(commandToExecute);
        else
            h = COMSERVER_HANDLE;
            Execute(h, commandToExecute);
            d = GetWorkspaceData(h, 'd', 'base');
        end
        % End for using java

    catch
        % In the worst case return only C:\
        errorMessage = lasterr;
        d = {'C:\'};
    end

else
    d = {filesep};
end

% Set FileRootsData
setappdata(0, 'FileRootsData', d);