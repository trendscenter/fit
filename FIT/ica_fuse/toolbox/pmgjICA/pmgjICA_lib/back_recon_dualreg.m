function varargout=back_recon_dualreg(sesInfo, dataMat)
    
    if isfield(sesInfo, 'outputDir')
            outputDir = sesInfo.outputDir;
        else
            outputDir=sesInfo.userInput.pwd;
    end

    % Get modality, data title and compset fields
    [modalityType, dataTitle, compSetFields] = mygift_get_modality;
    compSetFields = {'ic_gmcm','ic_gm', 'ic_icn', 'tc_gmcm','tc_gm', 'tc_icn'};

    icaStr = mygift_icaAlgorithm;
    algorithmName = deblank(icaStr(sesInfo.algorithm, :));

    icain =[sesInfo.ica_mat_file, '.mat'];  %%%% W from tst_ica.mat
    if isfield(sesInfo, 'backReconType')
            backReconType = sesInfo.backReconType;
    end

    backReconOptions = cellstr(mygift_backReconOptions);
    backReconInd = strmatch(backReconType, lower(backReconOptions), 'exact');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load(fullfile(outputDir, icain), 'icasig');  %%% icasig: 2 143733
    numComp=size(icasig, 1);

    preproc_type = sesInfo.preproc_type;
    backReconType = sesInfo.backReconType;
    back_reconstruction_mat_file = sesInfo.back_reconstruction_mat_file;
    outputDir = sesInfo.outputDir;
    diffTimePoints = sesInfo.diffTimePoints;
    numOfSub = sesInfo.numOfSub;
    numOfSess = sesInfo.numOfSess;
    num_fmri=sesInfo.userInput.numOffMRI;
    verbose = 0;

    icasig = icasig';
    dataGM=dataMat(1:length(sesInfo.userInput.maskmod1), 1:num_fmri:size(dataMat, 2));
    dataICN=dataMat(length(sesInfo.userInput.maskmod1)+1:end, :);
    
    meanMap_gmcm = zeros(numComp, length(sesInfo.userInput.mask_joint));
    meanMap_gm = zeros(numComp, length(sesInfo.userInput.maskmod1));
    meanMap_cm = zeros(numComp, length(sesInfo.userInput.maskmod2));
    
    
    for n_icn=1:num_fmri
        

        s_icn=dataICN(:, n_icn:num_fmri:size(dataICN, 2));  %%% ICN1, 2, 3, .... , 53
        sgcomsig=[];
        
        for n_sub=1:numOfSub
            compsig=[dataGM(:,n_sub); s_icn(:,n_sub)];
            %%size(compsig)
            sgcomsig=[sgcomsig,compsig];
            %%size(sgcomsig) %%% VoxelsxSubjects       
        end
        tmpCompSetFields = compSetFields;
        %%% for seperately
        dataGm1=sgcomsig(1:length(sesInfo.userInput.maskmod1), :);   
        dataCm1=sgcomsig(length(sesInfo.userInput.maskmod1)+1:end, :);  %%%whitesigCm';
         % Call pre-processing function
        if ~strcmpi(preproc_type, 'none')
            sgcomsig = mygift_preproc_data(sgcomsig, preproc_type, verbose);
            dataGm1 = mygift_preproc_data(dataGm1, preproc_type, verbose);
            dataCm1 = mygift_preproc_data(dataCm1, preproc_type, verbose);
        end

         %% Dual regression
         %% tc:  subxcomp
         %% spatial_maps:  compxvoxels

         [tc_gmcm, spatial_maps_gmcm] = mygift_dual_regress(sgcomsig, icasig);
         [tc_gm, spatial_maps_gm] = mygift_dual_regress(dataGm1, icasig(1:length(sesInfo.userInput.maskmod1), :));
         [tc_cm, spatial_maps_cm] = mygift_dual_regress(dataCm1, icasig(length(sesInfo.userInput.maskmod1)+1:end, :));

         %%disp('Dual regression output')
         %%size(tc_gmcm)  %%%% 4 2
         %%size(spatial_maps_gmcm)   %%%%2x143733

         compSet = struct(tmpCompSetFields{1}, spatial_maps_gmcm, tmpCompSetFields{2}, spatial_maps_gm, tmpCompSetFields{3}, spatial_maps_cm,...
             tmpCompSetFields{4}, tc_gmcm,  tmpCompSetFields{5}, tc_gm, tmpCompSetFields{6}, tc_cm);
         %%%compSet_gm = struct(tmpCompSetFields{1}, spatial_maps_gm, tmpCompSetFields{2}, tc_gm);
         %%%compSet_cm = struct(tmpCompSetFields{1}, spatial_maps_cm, tmpCompSetFields{2}, tc_cm);

         meanMap_gmcm = meanMap_gmcm + spatial_maps_gmcm;
         meanMap_gm = meanMap_gm + spatial_maps_gm;
         meanMap_cm = meanMap_cm + spatial_maps_cm;

          %%% Save Results          
          subFile = [back_reconstruction_mat_file, num2str(n_icn), '.mat'];
          txtMsg = ['-saving back reconstructed ica data for set ', num2str(n_icn), ' -> ',subFile];
          disp(txtMsg);
          mygift_parSave(fullfile(outputDir, subFile), {compSet}, {'compSet'});
          %%mygift_parSave(fullfile(outputDir, subFile), {compSet_gm}, {'compSet_gm'});
          %%mygift_parSave(fullfile(outputDir, subFile), {compSet_cm}, {'compSet_cm'});

         %%if n_icn>2
         %%    break
         %%end 

    end

    meanMap_gmcm = meanMap_gmcm./num_fmri;
    meanMap_gm = meanMap_gm./num_fmri;
    meanMap_cm = meanMap_cm./num_fmri;

    meanfile=fullfile(outputDir, 'meanmaps.mat');
    mygift_save(meanfile, 'meanMap_gmcm', 'meanMap_gm', 'meanMap_cm');

    if (nargout > 0)
        varargout{1} = meanMap_gmcm;
        varargout{2} = meanMap_gm;
        varargout{3} = meanMap_cm;
        varargout{4} = compSet_gmcm;
        varargout{5} = compSet_gm;
        varargout{6} = compSet_cm;
    else
        varargout = {};
    end







    

