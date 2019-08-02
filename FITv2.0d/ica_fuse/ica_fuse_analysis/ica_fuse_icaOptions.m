function [ICA_Options] = ica_fuse_icaOptions(dataSize, algorithm_index, handle_visibility)
% Return arguments depending on the ICA algorithm
% ICA Options is a cell arary with strings followed by values

ICA_Options = {};

if ~exist('handle_visibility', 'var')
    handle_visibility = 'on';
end

% determine the data size
numComps = dataSize(1);
frames = dataSize(2);
%[numComps, frames] = size(data);

icaAlgo = ica_fuse_icaAlgorithm; % available ICA algorithms

if ischar(algorithm_index)
    algorithm_index = strmatch(lower(algorithm_index), lower(icaAlgo), 'exact');
    if isempty(algorithm_index)
        error('Check the algorithm name you specified');
    end
end


% check the number of options
if algorithm_index <= size(icaAlgo, 1)
    % selected ica algorithm
    ica_algorithm = deblank(icaAlgo(algorithm_index , :));
else
    disp(['Presently there are ', num2str(size(icaAlgo, 1)), ' algorithms']);
    ICA_Options = {};
    return;
end

DEFAULT_BIAS = 'on';

global FLAG_STANDARDIZE_SUB;

if ~isempty(FLAG_STANDARDIZE_SUB)
    if strcmpi(FLAG_STANDARDIZE_SUB, 'yes')
        DEFAULT_BIAS = 'off';
    end
end

% Return options for the Infomax algorithm or Optimal ICA algorithm
if strcmpi(ica_algorithm, 'infomax') | strcmpi(ica_algorithm, 'sdd ica')
    
    % Defaults for the inofmax algorithm
    MAX_WEIGHT = 1e8;       % guess that weights larger than this have blown up
    DEFAULT_BLOCK = floor(sqrt(frames/3));
    DEFAULT_weight = 0.0;
    DEFAULT_LRATE = 0.015/log(numComps);
    DEFAULT_ANNEALSTEP   = 0.9;
    DEFAULT_ANNEALDEG = 60;
    DEFAULT_STOP = 1e-6;
    DEFAULT_MAXSTEPS = 512;
    DEFAULT_MOMENTUM = 0.0;
    DEFAULT_EXTENDED = 0.0;
    
    % dialog Title
    dlg_title = 'Select the Options for the Infomax algorithm';
    
    numParameters = 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = [['Select block less than ', num2str(frames)]  ...
        [' where Default = ', num2str(DEFAULT_BLOCK)]];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_BLOCK);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'block';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = ['Select stop where Default =', num2str(DEFAULT_STOP)];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_STOP);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'stop';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = ['Select weight where Maxweight =', num2str(MAX_WEIGHT)];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_weight);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'weights';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString ='Select lrate where min = 0.000001 and max = 0.1';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_LRATE);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'lrate';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;  % get the characters associated with the infomax algorithm
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = ['Select maxsteps where Default =' num2str(DEFAULT_MAXSTEPS)];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_MAXSTEPS);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'maxsteps';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Select anneal between (0 1]';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_ANNEALSTEP);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'anneal';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Select annealdeg between [0 180]';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_ANNEALDEG);;
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'annealdeg';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select momentum between [0 1]';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_MOMENTUM);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'momentum';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = ['Select Extended where default = ' num2str(DEFAULT_EXTENDED)];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(DEFAULT_EXTENDED);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'extended';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    %     inputText(numParameters).promptString = ['Select PCA in the range [1:', num2str(numComps - 1), ']'];
    %     inputText(numParameters).uiType = 'edit';
    %     inputText(numParameters).answerString = num2str(numComps - 1);
    %     inputText(numParameters).dataType = 'numeric';
    %     inputText(numParameters).tag = 'pca';
    %
    %     numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Number of Components';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = num2str(numComps);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'ncomps';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = ' Select Posact';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'on', 'off'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'posact';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = ' Select Sphering';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'on', 'off'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'sphering';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 2;
    
    
    numParameters = numParameters + 1;
    
    biasString = {'on', 'off'};
    matchIndex = strmatch(lower(DEFAULT_BIAS), lower(biasString), 'exact');
    
    
    inputText(numParameters).promptString = 'Select Bias';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = biasString;
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'bias';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = matchIndex;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select Verbose';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'on', 'off'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'verbose';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    
    % Select options for the Fast ICA algorithm
elseif strcmpi(ica_algorithm, 'fast ica')
    
    % Defaults for the FastICA
    epsilon = 0.0001;
    %epsilon = 1e-6;
    maxNumIterations  = 1000;
    maxFinetune       = 5;
    sampleSize        = 1;
    
    % dialog Title
    dlg_title = 'Select the Options for the Fast ICA algorithm';
    
    numParameters = 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = ['Select epsilon where default = ', num2str(epsilon)];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(epsilon);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'epsilon';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString =  ['Select maximum iterations where default = ', num2str(maxNumIterations)];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(maxNumIterations);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'maxNumIterations';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = ['Select maximum finetune where default = ', num2str(maxFinetune)];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(maxFinetune);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'maxFinetune';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = ['Select sampleSize between [0, 1] where default = ', num2str(sampleSize)];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(sampleSize);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'sampleSize';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;  % get the characters associated with the infomax algorithm
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Number of components';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = num2str(numComps);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'numOfIC';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Select verbose as on/off';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'on', 'off'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'verbose';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Select interactivePCA as on/off';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'off', 'on'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'interactivePCA';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select approach as defl/symm';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'defl', 'symm'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'approach';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 2;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select finetune from one of the following:  [off, pow3, tanh, gauss, skew]';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'off', 'pow3', 'tanh', 'gauss', 'skew'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'finetune';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select stabilization as on/off';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'off', 'on'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'stabilization';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select displayMode as on/off';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'off', 'on'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'displayMode';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select g from one of the following: [pow3, tanh, gauss, skew]';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'pow3', 'tanh', 'gauss', 'skew'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'g';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 2;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select Only from the following:';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'all', 'pca', 'white'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'only';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
elseif strcmpi(ica_algorithm, 'ccica')
    % Select options for the CCICA algorithm
    
    %% Options for CCICA
    % Maximum no. of constraint steps
    MAX_CONSTRAINT_STEPS = 10;
    % Min correlation between un-mixing matrix before and after adding the
    % constraint
    if frames>10000
        MIN_CORRELATION = 0.95;
    else
        MIN_CORRELATION=0.98;
    end
    %the number of components to be contrained
    NUM_CONSTRAINED_COMP = numComps;
    
    % dialog Title
    dlg_title = 'Select the options for the CCICA algorithm';
    
    numParameters = 1;
    
    % Max no. of constraint steps
    inputText(numParameters).promptString = 'Select the no. of constraint steps in the range [8-15]';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(MAX_CONSTRAINT_STEPS);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'max_constraint_steps';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % Min Correlation
    inputText(numParameters).promptString = 'Select the min correlation in the range [0.95-0.99] between un-mixing matrix before and after adding constraint';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(MIN_CORRELATION);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'min_correlation';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    % No. of constrained components
    inputText(numParameters).promptString = ['Select the no. of constrained components in the range [', ...
        num2str(ceil(NUM_CONSTRAINED_COMP/2)), '-', num2str(NUM_CONSTRAINED_COMP), ']'];
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(NUM_CONSTRAINED_COMP);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'num_constrained_components';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    
elseif strcmpi(ica_algorithm, 'erbm')
    %% ERBM
    
    dlg_title = 'Select the options for the ERBM algorithm';
    
    numParameters = 1;
    
    inputText(numParameters).promptString = 'Enter filter length';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(min(11, frames/50));
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'filter_length';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
elseif strcmpi(ica_algorithm, 'iva-g')
    %% IVA-G
    
    dlg_title = 'Select the options for the IVA-G algorithm';
    
    numParameters = 1;
    
    inputText(numParameters).promptString = 'Select approach for second order IVA';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = char('newton', 'quasi');
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'opt_approach';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Enter max no of iterations';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(512);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'maxIter';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Enter stopping tolerance';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '1e-6';
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'WDiffStop';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Enter initial step size scaling';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '1.0';
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'alpha0';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select verbose';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = char('Yes', 'No');
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'verbose';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
elseif strcmpi(ica_algorithm, 'iva-ggd')
    %% IVA-GGD
    
    dlg_title = 'Select the options for the IVA-GGD algorithm';
    
    numParameters = 1;
    
    inputText(numParameters).promptString = 'Enter max no of iterations';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(512);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'maxIter';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Enter stopping tolerance';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '1e-6';
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'termThreshold';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Enter initial step size scaling';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '1.0';
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'alpha0';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 0;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select verbose';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = char('Yes', 'No');
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'verbose';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
end

if exist('inputText', 'var')
    
    % Input dialog box
    answer = ica_fuse_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', handle_visibility);
    
    % ICA options with flags and the values corresponding to it
    ICA_Options = cell(1, 2*length(answer));
    
    for i = 1:length(answer)
        ICA_Options{2*i - 1} = inputText(i).tag;
        ICA_Options{2*i} = answer{i};
    end
    
end