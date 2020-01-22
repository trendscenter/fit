function [compData, anatData, meanData, groupCompData, meanDataLegend, groupCompLegend, ...
    text_left_right, plotType, HInfo] = ica_fuse_loadCompData(varargin)
% load component data


ica_fuse_defaults;
global DETREND_NUMBER;
global FLIP_ANALYZE_IM;

interpMessage = 'Interpolating components. Please wait ...';
interpTitle = 'Interpolating Components';

oldDir = pwd;
meanData = [];
b = [];
HInfo = [];
meanDataLegend = {};
groupCompLegend = {};
groupCompData = [];
groupCompFiles = [];
flip_analyze_images = [];
modalityName = [];
% Loop over number of arguments
for ii = 1:2:nargin
    % Get the required vars
    if strcmpi(varargin{ii}, 'component_files')
        % component files
        component_files = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'anatomical_file')
        % anatomical file
        anatomical_file = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'component_numbers')
        % component numbers
        component_numbers = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'anatomical_view')
        % anatomical view
        anatomicalView = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'slices_in_mm')
        % slices in mm
        slices_in_m = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'interp_message')
        interpMessage = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'interp_title')
        interpTitle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'input_files')
        inputFiles = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'mask_ind')
        mask_ind = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'feature_info')
        featureInfo = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'voxels')
        voxels = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'outputdir')
        outputDir = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'group_component_files')
        groupCompFiles = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'flip_analyze_images')
        flip_analyze_images = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'modality')
        modalityName = varargin{ii + 1};
    end
    % end for getting the required vars
    
end
% end loop over number of arguments

text_left_right = [];
anatData = [];

% Get the number of components
[countTimePoints] = ica_fuse_get_countTimePoints(component_files);

maxComp = max(component_numbers);

if maxComp > countTimePoints
    error(['Component number requested exceeds the number of actual components']);
end

% First component file
firstFile = deblank(component_files(1, :));

% Get the extension and the file path
[compPath, firstFname, extn] = fileparts(firstFile);

% Load the component data appropriately
if strcmpi(extn, '.img') || strcmpi(extn, '.nii')
    
    %%%%%%%% Display warning message %%%%%%%%%%
    compV = ica_fuse_getVol(firstFile, 1);
    classComp = 'analyze';
    if strcmpi(class(compV(1).private), 'ica_fuse_nifti')
        classComp = 'nifti';
    end
    
    [status] = checkMatFile(firstFile);
    
    if ~isempty(flip_analyze_images)
        if strcmpi(classComp, 'analyze') && (status == 0)
            if (flip_analyze_images ~= FLIP_ANALYZE_IM)
                warning(['Flip parameter (', num2str(FLIP_ANALYZE_IM), ...
                    ') specified in ica_fuse_defaults is different from the one specified during the analysis (', ...
                    num2str(flip_analyze_images), ')']);
            end
        end
    end
    %%%%%%%% End for displaying warning message %%%%%%%%%%
    
    plotType = 'image';
    try
        helpHandle = helpdlg(interpMessage, interpTitle);
    catch
    end
    % Get the structural volume
    structVol = ica_fuse_get_vol_nifti(anatomical_file);
    [images, coords, HInfo, slices, text_left_right] = ica_fuse_resizeImage(structVol, component_files, anatomicalView, ...
        slices_in_m, component_numbers);
    % Get the component data
    compData = images(2:end, :, :, :);
    anatData = images(1, :, :, :);
    
    anatData = squeeze(anatData);
    if length(size(anatData)) == 2
        anatData = reshape(anatData, [size(anatData), 1]);
    end
    
    
    try
        delete(helpHandle);
    catch
    end
    
else
    
    
    
    plotType = 'timecourse';
    
    if (~isempty(modalityName) && strcmpi(modalityName, 'gene'))
        plotType = 'SNP';
    end
    
    % Treat the file as ascii
    [compData] = ica_fuse_loadData(component_files(component_numbers, :));
    
    if (strcmpi(plotType, 'timecourse'))
        
        % Number of groups
        numGroups = size(featureInfo.groupNames, 1);
        
        % Number of subjects
        numSubjects = featureInfo.numSubjects;
        
        % Input data
        input_data = ica_fuse_loadData(inputFiles(1:numSubjects(1), :));
        
        % Mean of Y data
        meanData = mean(squeeze(input_data(mask_ind, 2, :)), 2);
        
        meanDataLegend{1} = ['Mean of ', deblank(featureInfo.groupNames(1, :))];
        
        clear input_data;
        
        % Interpolation factor
        interpFactor = ceil(voxels / size(meanData, 1));
        
        if interpFactor ~= 0
            % Resample mean data
            meanData = ica_fuse_resample(meanData, voxels, size(meanData, 1));
            
            %%%%%%%%%%%% Resample component data %%%%%%%%%%%%%%%%%%
            [compData] = ica_fuse_resampleCompData(compData, length(meanData), voxels, length(mask_ind));
            
        end
        
        % Replicate meanData over groups
        meanData = repmat(meanData, 1, numGroups);
        
        if ~isempty(groupCompFiles)
            groupCompLegend{1} = ['Group ', deblank(featureInfo.groupNames(1, :))];
            groupCompData = zeros([size(compData, 1), size(compData, 2), length(component_numbers), numGroups]);
            if length(component_numbers) > 1 & numGroups > 1
                %groupCompData(:, :, :, 1) = ica_fuse_loadData(str2mat(groupCompFiles(1).comp.name), component_numbers);
                tempCompData = ica_fuse_loadData(str2mat(groupCompFiles(1).comp.name), component_numbers);
                groupCompData(:, :, :, 1) = ica_fuse_resampleCompData(tempCompData, size(meanData, 1), voxels, ...
                    length(mask_ind));
                clear tempCompData;
                
            else
                groupCompData = ica_fuse_loadData(str2mat(groupCompFiles(1).comp.name), component_numbers);
                groupCompData = ica_fuse_resampleCompData(groupCompData, size(meanData, 1), voxels, length(mask_ind));
            end
        end
        
        % check if the number of groups is greater than 1
        if numGroups > 1
            startInd = numSubjects(1) + 1;
            % Loop over groups
            for nn = 2:numGroups
                endInd = sum(numSubjects(1:nn));
                input_data = ica_fuse_loadData(inputFiles(startInd:endInd, :));
                tempData = mean(squeeze(input_data(mask_ind, 2, :)), 2);
                clear input_data;
                if interpFactor ~= 0
                    meanData(:, nn) = ica_fuse_resample(tempData, voxels, size(tempData, 1));
                end
                clear tempData;
                meanDataLegend{nn} = ['Mean of ', deblank(featureInfo.groupNames(nn, :))];
                if ~isempty(groupCompFiles)
                    tempCompData = ica_fuse_loadData(str2mat(groupCompFiles(nn).comp.name), component_numbers);
                    groupCompData(:, :, :, nn) = ica_fuse_resampleCompData(tempCompData, size(meanData, 1), voxels, ...
                        length(mask_ind));
                    clear tempCompData;
                    %groupCompData(:, :, :, nn) = ica_fuse_loadData(str2mat(groupCompFiles(nn).comp.name), ...
                    %    component_numbers);
                    groupCompLegend{nn} = ['Group ', deblank(featureInfo.groupNames(nn, :))];
                end
                startInd = endInd + 1;
            end
            % End loop over groups
        end
        % End for checking if the number of groups is greater than 1
        
        % Get the data for the selected component numbers
        %compData = compData(:, :, component_numbers);
        % Make component number as the first dimension
        compData = permute(compData, [3, 1, 2]);
        
    else
        compData = permute(compData, [3, 1, 2]);
        compData = compData(:, :, 2);
        
    end
    
end

function [status] = checkMatFile(file_name)
% Check MAT file name and return 1 if it exists

status = 0;
[pp, fName, extn] = fileparts(file_name);

% MAT file name
if strcmp(extn, '.img')
    matFileName = fullfile(pp, [fName, '.mat']);
else
    matFileName = fullfile(pp, [fName, '.MAT']);
end

status = exist(matFileName, 'file');

if status > 1
    status = 1;
end