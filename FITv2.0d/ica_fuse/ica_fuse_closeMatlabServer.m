function ica_fuse_closeMatlabServer
% close Matlab server on PC's for nojvm

global COMSERVER_HANDLE;

try
    % Quit COM server only for windows
    if ispc
        if ~usejava('jvm')
            disp('Closing Matlab Server ...');
            fprintf('\n');
            Quit(COMSERVER_HANDLE);
            delete(COMSERVER_HANDLE);
            clear global COMSERVER_HANDLE;
        end
    end

catch

end