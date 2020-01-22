function fusionInfo = ica_fuse_dataReduction(fusionInfo, data)
%% Data reduction step for fusion analysis. This information will be stored
% in PCA MAT files. The MAT file name is dependent on the global variable
% DATA_REDUCTION_FILE (see ica_fuse_defaults) and combination number.
%
% Inputs:
% 1. fusionInfo - Fusion Info data structure
% 2. data - Stacked data
%
% Outputs:
% fusionInfo - Update Fusion Info data structure
%

% Load defaults`
ica_fuse_defaults;
global DATA_REDUCTION_FILE;


%% Output directory
outputDir = fusionInfo.run_analysis.outputDir;
combinationName = fusionInfo.run_analysis.currentCombName;
comb_number = fusionInfo.run_analysis.currentComb;
numComp = fusionInfo.run_analysis.numComp;
allComb = fusionInfo.run_analysis.all_comb;

%% Get type of PCA and reference vector
type_pca = 'standard';
reference = [];
try
    type_pca = fusionInfo.run_analysis.type_pca;
    reference = fusionInfo.run_analysis.reference;
catch
    
end

disp('-----------------------------------------------------------------------------------------------');
disp(['Doing data reduction for ',  combinationName]);
disp('-----------------------------------------------------------------------------------------------');

fprintf('\n');

%% Save data reduction information in file
pcaFile = [fusionInfo.run_analysis.prefix, DATA_REDUCTION_FILE, 'comb_', num2str(comb_number), '.mat'];
fusionInfo.run_analysis.pcaFiles(comb_number).name = pcaFile;
fusionInfo.run_analysis.pcaFiles(comb_number).combinationName = combinationName;
pcaFile = fullfile(outputDir, pcaFile);

featureNames = cellstr(char(fusionInfo.run_analysis.dataInfo(1).feature.name));

ica_algo = ica_fuse_icaAlgorithm;
ica_algo = cellstr(ica_algo);
ica_type = ica_algo{fusionInfo.run_analysis.algorithm};
modalities = cellstr(char(fusionInfo.run_analysis.dataInfo(1).feature.modality));
modalities = modalities(allComb{comb_number});
dims = fusionInfo.run_analysis.newDims(allComb{comb_number});

if (strcmpi(ica_type, 'iva-g') || strcmpi(ica_type, 'iva-ggd'))
    
    V = cell(1, length(modalities));
    Lambda = V;
    whitesig = V;
    whiteM = V;
    dewhiteM = V;
    endM = 0;
    for nM = 1:length(modalities)
        %% Do PCA and whitening
        startM = endM + 1;
        endM = endM + dims(nM);
        dat  = data(:, startM:endM)';
        if (~strcmpi(modalities{nM}, 'gene'))
            dat = ica_fuse_remove_mean(dat);
        end
        [V{nM}, Lambda{nM}, whitesig{nM}, whiteM{nM}, dewhiteM{nM}] = ica_fuse_calculate_pca(dat, numComp);
    end
    
    ica_fuse_save(pcaFile, 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM', 'combinationName');
    
else
    
    
    
    if (strcmpi(type_pca, 'standard') || strcmpi(type_pca, 'reference'))
        
        %% Do PCA and whitening
        [V, Lambda, whitesig, whiteM, dewhiteM] = ica_fuse_calculate_pca(data, numComp, type_pca, reference);
        
        ica_fuse_save(pcaFile, 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM', 'combinationName');
        
        
    else
        
        if (~strcmpi(type_pca, 'mccar'))
            try
                numPC = fusionInfo.run_analysis.cca_opts.numPC;
            catch
                numPC = [fusionInfo.run_analysis.cca_opts.numPC1, fusionInfo.run_analysis.cca_opts.numPC2];
            end
        end
        
        if (strcmpi(type_pca, 'cca') && all(fusionInfo.run_analysis.newDims == fusionInfo.run_analysis.newDims(1)))
            
            
            disp('Doing individual PCAs on the features ...');
            
            %% Do PCA on feature 1
            [V1, Lambda1, whitesig1, whiteM1] = preCCA(data(:, 1:fusionInfo.run_analysis.newDims(1)), numPC(1));
            
            %% Do PCA 2 on feature 2
            [V2, Lambda2, whitesig2, whiteM2] = preCCA(data(:, fusionInfo.run_analysis.newDims(1) + 1:end), numPC(2));
            
            disp('Doing CCA on the PCA reduced matrices ...');
            
            %% PCA on the correlation matrices of feature 1 and 2
            [V, Lambda, whitesig, whiteM] = ica_fuse_calculate_pca(whitesig1*whitesig2', numComp);
            D = sum(whitesig(1, :).^2);
            E1 = V';
            E2 = whiteM*(whitesig1*whitesig2')/sqrt(D);
            Lambda = sqrt(Lambda*D);
            whitesig = [E1*whitesig1, E2*whitesig2];
            whiteM = {E1*whiteM1, E2*whiteM2};
            dewhiteM = {pinv(whiteM{1}), pinv(whiteM{2})};
            
            ica_fuse_save(pcaFile, 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM', 'E1', 'E2', 'combinationName');
            
        elseif (strcmpi(type_pca, 'mccar'))
            %% MCCAR
            %opts.numc = numComp;
            disp('Doing MCCAR on the features ...');
            endM = 0;
            dat = cell(1, length(modalities));
            for nM = 1:length(modalities)
                startM = endM + 1;
                endM = endM + dims(nM);
                dat{nM}  = data(:, startM:endM)';
                if (~strcmpi(modalities{nM}, 'gene'))
                    dat{nM}  = ica_fuse_remove_mean(dat{nM});
                end
            end
            
            [dewhiteM, whitesig, whiteM] = ica_fuse_mcca_reference(dat, numComp, reference);
            ica_fuse_save(pcaFile, 'whitesig', 'dewhiteM', 'whiteM', 'combinationName');
            
        else
            
            %% Use MCCA
            disp('Doing MCCA on the features ...');
            
            if (fusionInfo.run_analysis.numFeatures == 2)
                [C1, C2, E1, E2, Lambda] = ica_fuse_mcca(data(:, 1:fusionInfo.run_analysis.newDims(1)), data(:, fusionInfo.run_analysis.newDims(1) + 1:end), ...
                    numPC(1), numPC(2), numComp, cellstr(char(fusionInfo.run_analysis.dataInfo(1).feature.name)));
                
                whitesig = [C1, C2];
                dewhiteM = {data(:, 1:fusionInfo.run_analysis.newDims(1))*pinv(C1), data(:, fusionInfo.run_analysis.newDims(1) + 1:end)*pinv(C2)};
                whiteM = {pinv(dewhiteM{1}), pinv(dewhiteM{2})};
                ica_fuse_save(pcaFile, 'Lambda', 'whitesig', 'whiteM', 'dewhiteM', 'combinationName');
            else
                [whitesig, dewhiteM, whiteM] = ica_fuse_mccaN(data, fusionInfo.run_analysis.newDims, numPC, numComp, featureNames);
                ica_fuse_save(pcaFile, 'whitesig', 'whiteM', 'dewhiteM', 'combinationName');
            end
            
        end
        
        
    end
    
end

disp(['Data reduction information for ', combinationName, ' is saved in ', pcaFile]);

fprintf('\n');

%% Save fusion file
fusionFile = fusionInfo.run_analysis.fusionFile;
ica_fuse_save(fusionFile, 'fusionInfo');

disp('-----------------------------------------------------------------------------------------------');
disp(['Done data reduction for ',  combinationName]);
disp('-----------------------------------------------------------------------------------------------');

fprintf('\n');


function [V, Lambda, whitesig, whiteM] = preCCA(dat, numComp)
%% Do PCA and normalize data before doing CCA
%

%dat = ica_fuse_remove_mean(dat')';
[V, Lambda, whitesig, whiteM] = ica_fuse_calculate_pca(dat, numComp);
D = sum(whitesig(1, :).^2);
whitesig = whitesig / sqrt(D);
whiteM = whiteM / sqrt(D);