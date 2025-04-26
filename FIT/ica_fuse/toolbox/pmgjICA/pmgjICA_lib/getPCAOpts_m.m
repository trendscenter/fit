function pca_opts = getPCAOpts_m(sesInfo)
    % Get PCA opts

    
    
    isOptsCell = 0;
    try
        isOptsCell = iscell(sesInfo.pca_opts);
    catch
    end
    
    if (~isOptsCell)
        
        if (isfield(sesInfo, 'pca_opts'))
            tmp_pca_opts = sesInfo.pca_opts;
        end
        
        %%% Group PCA Options
        pcaType = 'standard';
        if (isfield(sesInfo, 'pcaType'))
            pcaType = sesInfo.pcaType;
        end
        
        sesInfo.pcaType = pcaType;
        sesInfo = mygift_check_pca_opts(sesInfo);
        pcaType = sesInfo.pcaType;
        tmp_pca_opts = sesInfo.pca_opts;
        
        pca_opts{1} = struct('pcaType', pcaType, 'pca_opts', tmp_pca_opts);
        
    else
        
        pca_opts = sesInfo.pca_opts;
        for nP = 1:length(pca_opts)
            pca_opts{nP} = mygift_check_pca_opts(pca_opts{nP});
        end
        
    end
    
    if (length(pca_opts) ~= sesInfo.numReductionSteps)
        pca_opts = repmat(pca_opts(1), 1, sesInfo.numReductionSteps);
    end