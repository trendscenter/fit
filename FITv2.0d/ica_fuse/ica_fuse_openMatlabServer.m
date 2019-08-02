function ica_fuse_openMatlabServer
% open Matlab server on PC's for nojvm

global COMSERVER_HANDLE;

if ispc
    % no java
    if ~usejava('jvm')
        % check Matlab server handle
        if ~(checkCOM(COMSERVER_HANDLE))
            % open matlab COM automation server
            comOpenStr = 'Opening Matlab server ...';
            disp(comOpenStr);
            COMSERVER_HANDLE = actxserver('matlab.application');
            set(COMSERVER_HANDLE, 'visible', 0); % make it not visible          
        end
        % check if it is Matlab server handle

    end
    % no java
end
% for pc

function status = checkCOM(h)
% Detect if it is COM object

status = 0;

if ~isempty(h)
    % MATLAB version
    versionNum = ica_fuse_get_matlab_version;

    if (versionNum < 14)
        status = isa(h, 'COM.matlab.application');
    else
        status = iscom(h);
    end
end