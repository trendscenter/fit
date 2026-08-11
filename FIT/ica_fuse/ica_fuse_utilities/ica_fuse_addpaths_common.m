function ica_fuse_addpaths_common()
    %icatb_ADDPATHS_COMMON Summary of this function goes here
    b_add_common = false;
    try
        if isempty(which('trd_util_slash.m'))
            b_add_common = true;
        end

        if b_add_common
            s_up1 = fileparts(which('fusion.m'));
            s_up2 = fileparts(s_up1);
            s_root = fileparts(s_up2);
        
            s_prefix = [s_root filesep 'code/common'];
            allDirs = strsplit(path, pathsep);
            if ~any(strncmp(allDirs, s_prefix, length(s_prefix)))
                addpath(genpath(s_prefix), '-end');
            end

        end

    catch
        disp('Warning: Common dir may be missing.');
    end