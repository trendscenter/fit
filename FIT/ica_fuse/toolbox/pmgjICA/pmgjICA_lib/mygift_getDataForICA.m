function data = mygift_getDataForICA(sesInfo, algorithmName, dataMat, reduction_sep, statusHandle)
%% Load data before applying ICA (reduced data from pca excluding spatially constrained ica algorithms)
%


if ~exist('nReduc', 'var')
    nReduc=1;
end

if (ischar(sesInfo))
    load(sesInfo);
end


if (~exist('statusHandle', 'var'))
    statusHandle = [];
end

modalityType = 'fMRI';

try

    if isfield(sesInfo, 'modality')
        modalityType = sesInfo.modality;
    else
        modalityType = sesInfo.userInput.modality;
    end
    
catch

end

if isfield(sesInfo, 'outputDir')
    outputDir = sesInfo.outputDir;
else
    outputDir = sesInfo.userInput.pwd;
end

mask_ind = sesInfo.mask_ind;

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    try
        useTemporalICA = strcmpi(sesInfo.group_ica_type, 'temporal');
    catch
    end
end

if (isempty(ica_fuse_findstr(lower(algorithmName),'iva')))
    if strcmpi(sesInfo.userInput.pjICA, 'yes') && sesInfo.userInput.groupica_algorithm==1
        %%%% Stacking the output of the first stage
        data_com=[];
        for n_icn=1:sesInfo.userInput.numOffMRI
            pcicafile=fullfile(outputDir, [sesInfo.data_reduction_mat_file, '_ica_1-', num2str(n_icn), '.mat']);
            load(pcicafile, 'icasig');
            data_com=[data_com; icasig];    
            clear icasig;

            if n_icn>2
               break
            end
            
        end
        %%size(data_com)  %%%% subject(ICN)xVoxels

    elseif strcmpi(sesInfo.userInput.pjICA, 'yes')
        if reduction_sep
            %try
                %pcain=cell(1, sesInfo.userInput.num_mod);
                %for nm=1:sesInfo.userInput.num_mod
                %    pcain{nm} = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(1), '_', num2str(nm),'.mat'];
                %end
            %catch
            disp('Signle file for multiple modalities..')
            pcain = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(1), '.mat'];  %%% 'tst_pca_r2-1.mat'
            %end
        else
            pcain = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(1), '_', num2str(nReduc),'.mat'];
        end
    else    
        pcain = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(1), '.mat'];  %%% 'tst_pca_r2-1.mat'
    end
    
    if (~useTemporalICA)

        if strcmpi(sesInfo.userInput.pjICA, 'yes')
            
            if sesInfo.userInput.groupica_algorithm==1
                data=data_com';
                clear data_com;

            elseif reduction_sep
                try
                    data=[];
                    for nm=1:sesInfo.userInput.num_mod
                        load(fullfile(outputDir, pcain{nm}));
                        data1=pcasig;
                        if size(data1, 1)== prod(sesInfo.userInput.HInfo_sep.DIM(1:3))
                            mask_indmod=sesInfo.userInput.(strcat('maskmod', num2str(nm)));
                            data1 = data1(mask_indmod, :);
                        end
                        data=[data; data1];  %%% VoxelsXSubjects
                        
                    end
                catch
                    disp('Read single file for multiple modalities..')                    
                    load(fullfile(outputDir, pcain));
                    data=pcasig;                

                end
            end

        else
            load(fullfile(outputDir, pcain));
            data = pcasig;
            if size(data, 1) == prod(sesInfo.HInfo.DIM(1:3))
                data = data(mask_ind, :);
            end
        end
    else
        if reduction_sep
            disp('Not yet implemented for the Temporal ICA');
            return
        else
            load(fullfile(outputDir, pcain), 'V');
            data = V;
            data = icatb_remove_mean(V);
        end
    end
    
    % Changed this code to allow the users add their own ICA algorithm
    % transpose data to equal components by volume
    data = data';
    
else
    if strcmpi(sesInfo.userInput.pjICA, 'yes')
        %%pcain = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(1), num2str(nReduc),'.mat'];
        fS = whos('-file', fullfile(outputDir, [sesInfo.data_reduction_mat_file, '1-1','_',num2str(nReduc), '.mat']));  %%% 'tst_pca_r1-1.mat'
    else
        fS = whos('-file', fullfile(outputDir, [sesInfo.data_reduction_mat_file, '1-1.mat']));  %%% 'tst_pca_r1-1.mat'
    end

    fNames = cellstr(char(fS.name));
    chkPcasig = isempty(strmatch('pcasig', fNames, 'exact'));
    
    disp('Stacking data across subjects ...');
    
    if (chkPcasig)
        
        if strcmpi(sesInfo.userInput.pjICA, 'yes')
            pcain=[sesInfo.data_reduction_mat_file, '1-1','_',num2str(nReduc), '.mat'];
        else
            pcain = [sesInfo.data_reduction_mat_file, '1-1.mat'];
        end
        info = load(fullfile(outputDir, pcain));
        VStacked = info.V;
        LambdaStacked = info.Lambda;
        
        for nD = 1:sesInfo.numOfSub*sesInfo.numOfSess
            if (nD == 1)
                startT = 1;
            else
                startT = sum(sesInfo.diffTimePoints(1:nD-1)) + 1;
            end
            endT = sum(sesInfo.diffTimePoints(1:nD));
            [whiteM, dewhiteM] = get_pca_info(VStacked(startT:endT, :), diag(LambdaStacked(nD, :)));

            if ~isempty(sesInfo.inputFiles)
                dat = icatb_read_data(sesInfo.inputFiles(nD).name, [], mask_ind);  %%%% 50485         220
            else
                dat = dataMat(:, ((nD-1)*num_fmri)+1:nD*num_fmri);
            end
            % Call pre-processing function
            if ~strcmpi(sesInfo.preproc_type, 'None')
                dat = icatb_preproc_data(dat, sesInfo.preproc_type, 0);
                % Remove mean per timepoint
                dat = icatb_remove_mean(dat, 0);
            end

            pcasig = dat*whiteM';
            
            if (nD == 1)
                if ~isempty(sesInfo.inputFiles)
                    data = zeros(sesInfo.numComp, length(mask_ind), sesInfo.numOfSub*sesInfo.numOfSess);
                else
                    data = zeros(sesInfo.numComp, size(dataMat, 1), sesInfo.numOfSub*sesInfo.numOfSess);
                end
            end
            
            data(:, :, nD) = pcasig';
            
            clear wM dat pcasig;
            
        end
        
    else
        
        % stack data of all subjects for doing IVA
        for nD = 1:sesInfo.numOfSub*sesInfo.numOfSess

            if strcmpi(sesInfo.userInput.pjICA, 'yes')
                pcain=[sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(nD), '_',num2str(nReduc), '.mat'];
            else
                pcain = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(nD), '.mat'];
            end

            load(fullfile(outputDir, pcain), 'pcasig');
            
            if (nD == 1)
                if ~isempty(sesInfo.inputFiles)
                    data = zeros(sesInfo.numComp, length(mask_ind), sesInfo.numOfSub*sesInfo.numOfSess);
                else
                    data = zeros(sesInfo.numComp, size(dataMat, 1), sesInfo.numOfSub*sesInfo.numOfSess);
                end
            end
            
            data(:, :, nD) = pcasig';
            
            clear pcasig;
            
            
        end
        
    end
    
end

function [whiteningMatrix, dewhiteningMatrix] = get_pca_info(V, Lambda)
%% Get Whitening and de-whitening matrix
%
% Inputs:
% 1. V - Eigen vectors
% 2. Lambda - Eigen values diagonal matrix
%
% Outputs:
% 1. whiteningMatrix - Whitening matrix
% 2. dewhiteningMatrix - Dewhitening matrix
%


whiteningMatrix = sqrtm(Lambda) \ V';
dewhiteningMatrix = V * sqrtm(Lambda);