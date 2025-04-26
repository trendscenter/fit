function data = loadData_m(files, mask_ind, preProcType, precisionType, varToLoad)
    %% Load data and preprocess if specified
    %

    if ~exist('mask_ind', 'var')
        mask_ind=[];
    end
    if ~exist('preProcType', 'var')
        preProcType='none';        
    end
    if ~exist('precisionType', 'var')
        precisionType='double';
    end
    if ~exist('varToLoad', 'var')
        varToLoad='pcasig';
    end
    
    
    data = cell(1, length(files));
    
    verbose = 0;
    if (length(files) > 1)
        verbose = 1;
    end
    
    for nF = 1:length(files)
        if (verbose)
            disp(['Loading data-set ', num2str(nF), ' ...']);
        end
        tmp = ica_fuse_read_data(files{nF}, [], mask_ind, precisionType, varToLoad);
        tmp = squeeze(tmp);
        size_d = size(tmp);
        tmp = reshape(tmp, prod(size_d(1:end-1)), size_d(end));
        if (~strcmpi(preProcType, 'none'))
            %Call pre-processing function
            tmp = ica_fuse_preproc_data(tmp, preProcType);
            if (~strcmpi(preProcType, 'remove mean per timepoint'))
                %Remove mean per timepoint
                tmp = ica_fuse_remove_mean(tmp, 1);
            end
        end
        data{nF} = tmp;
        
    end
    
    data = [data{:}];
