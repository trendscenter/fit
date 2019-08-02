function [answerString, answerValue] = ica_fuse_get_displayDefaults(controlName)

ica_fuse_defaults;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGE_VALUES;
global ANATOMICAL_PLANE;
global SORT_COMP;
global IMAGES_PER_FIGURE;
global ANATOMICAL_FILE;

controlName = lower(controlName);

switch controlName
    case 'image_values'
        answerString = str2mat('Positive and Negative', 'Positive', 'Absolute', 'Negative');
        [answerValue] = get_index(answerString, IMAGE_VALUES);
    case 'convert_to_z'
        answerString = str2mat('Yes', 'No');
        [answerValue] = get_index(answerString, CONVERT_TO_Z);
    case 'anatomical_plane'
        answerString = str2mat('Axial', 'Sagittal', 'Coronal');
        [answerValue] = get_index(answerString, ANATOMICAL_PLANE);
    case 'sort_comp'
        answerString = '';
        if ~isnumeric(SORT_COMP)
            answerValue = 0;
        else
            answerValue = SORT_COMP;
        end
    case 'z_threshold'
        answerString = num2str(Z_THRESHOLD);
        answerValue = 0;
    case 'images_per_figure'
        answerString = str2mat('1', '4', '9');
        [answerValue] = get_index(answerString, IMAGES_PER_FIGURE);
    otherwise
        answerString = '';
        answerValue = 0;
end


function [matchIndex] = get_index(promptString, index)

if ischar(index)
    
    matchIndex = strmatch(lower(index), lower(promptString), 'exact');
    if isempty(matchIndex)
        error('Check the input argument passed to image values');
    end
    
elseif isnumeric(index)
    
    if round(index) ~= index
        error('Index must be an integer');
    end
    
    if index < 0 | index > size(promptString, 1)
        error('Check the index parameter passed');
    end
    
    matchIndex = index;
    
else
    
    error('Index must be an integer or character');
    
end