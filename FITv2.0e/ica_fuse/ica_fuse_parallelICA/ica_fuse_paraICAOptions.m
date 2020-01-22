function [ICA_Options] = ica_fuse_paraICAOptions(numComp, handle_visibility, dims)
% ICA options for parallel ICA

CONSTRAINED_COMPONENTS = 3; % NUMBER OF COMPONENTS FROM EACH DATASET BEING CONSTRAINED
CONSTRAINED_CONNECTION = 0.3; % CORRELATION THRESHOLD TO BE CONSTRAINED; HIGH THRESHOLD WILL BE STRENGTHENED.
ENDURANCE = -1e-3; % the maximumlly allowed descending trend of entropy;

if numComp < CONSTRAINED_COMPONENTS
    CONSTRAINED_COMPONENTS = numComp;
end

if (length(numComp) > 1)
    numSubjects = numComp(2);
    numComp = numComp(1);
    CONSTRAINED_CONNECTION = p_to_r(numSubjects);
end

if ~exist('handle_visibility', 'var') || isempty(handle_visibility)
    handle_visibility = 'on';
end

useEndurance = 1;
if (exist('dims', 'var'))
    lrate = 0.015./log(dims);
    if (length(dims) == 3)
        useEndurance = 0;
    end
end

numParameters = 1;

% define all the input parameters in a structure
inputText(numParameters).promptString = 'Enter the maximum number of components from each data-set to be constrained.';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(CONSTRAINED_COMPONENTS);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'constrained_components';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 0;

numParameters = numParameters + 1;

% define all the input parameters in a structure
inputText(numParameters).promptString = 'Enter a correlation threshold to constrain components. Default is based on p < 0.05. Assume correction of such a value potentially meaningful.';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(CONSTRAINED_CONNECTION);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'constrained_connection';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 0;


if (useEndurance)
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Enter the maximumlly allowed descending trend of entropy.';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(ENDURANCE);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'endurance';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
else
    
    ANNEAL = 0.7;
    NOCHANGE = 1e-4;
    ANNEALCORR = 0.95;
    
    
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Select anneal.';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(ANNEAL, '%0.4f');
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'anneal';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Select correlation anneal step constant between (0, 1).';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(ANNEALCORR, '%0.4f');
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'annealcorr';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Select stopping tolerance.';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(NOCHANGE);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'nochange';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    
    
end


numParameters = numParameters + 1;

% define all the input parameters in a structure
inputText(numParameters).promptString = 'Enter max number of steps';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(512);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'maxsteps';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 0;

if (exist('lrate', 'var'))
    numParameters = numParameters + 1;
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Enter lrate for each modality';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(lrate);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'lrate';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
end


% Input dialog box
answer = ica_fuse_inputDialog('inputtext', inputText, 'Title', 'Select ICA options', 'handle_visibility', handle_visibility);

if isempty(answer)
    error('Parallel ICA options are not selected');
end

if answer{1} > numComp
    error('Error:numComp', 'Maximum number of components (%s) to be constrained exceeds the number of components (%s)', num2str(answer{1}), num2str(numComp));
end

if answer{2} > 1
    error('Correlation threshold to constrain components cannot be greater than 1');
end

% ICA options with flags and the values corresponding to it
ICA_Options = cell(1, 2*length(answer));

if ~isempty(answer)
    for i = 1:length(answer)
        ICA_Options{2*i - 1} = inputText(i).tag;
        ICA_Options{2*i} = answer{i};
    end
else
    disp('Using defaults options for parallel ICA ...');
    ICA_Options = {'constrained_components', CONSTRAINED_COMPONENTS, 'constrained_connection', CONSTRAINED_CONNECTION, 'endurance', ENDURANCE};
end


function r = p_to_r(N)
% p to r value

x = abs(ica_fuse_spm_invTcdf(0.05, N-2));
r = sqrt(1/((N-2) + x^2))*x;