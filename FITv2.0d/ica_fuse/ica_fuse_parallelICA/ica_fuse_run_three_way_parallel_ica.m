function [weights,sphere,activations,bias,signs,lrates,y,sim,idxt,num,sim_c,debug_entropy] = ica_fuse_run_three_way_parallel_ica(data,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)
%% %%%%%%%%%%%%%%%%%%%%% Sizes %%%%%%%%%%%%%%%%%%%%%%%%
% Sizes
data1=data{1};
data2=data{2};
data3=data{3};
%size and remember of data1
[chans(1),frames(1)]=size(data1);
urchans(1)=chans(1);
datalength(1)=frames(1);
%size and remember of data2
[chans(2),frames(2)]=size(data2);
urchans(2)=chans(2);
datalength(2)=frames(2);
%size and remember of data3
[chans(3),frames(3)]=size(data3);
urchans(3)=chans(3);
datalength(3)=frames(3);
%% %%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 1e-9;      % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      % anneal by multiplying lrate by this original 0,9 to 0.95 changed by JL
DEFAULT_ANNEALCORR   = 0.5;       % anneal by multiplying Crate by this original
DEFAULT_EXTANNEAL    = 0.98;      % or this if extended-ICA
DEFAULT_MAXSTEPS     = 512;       % ]top training after this many steps 512
DEFAULT_MOMENTUM     = 0.0;       % default momentum weight

DEFAULT_BLOWUP       = 1e9;       % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.9;       % if weights blowup, restart with lrate lower by this factor
MIN_LRATE            = 1e-8;       % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
DEFAULT_LRATE        = 0.015./log(chans);

% heuristic default - may need adjustment
%   for large or tiny data sets!
DEFAULT_BLOCK        = floor(sqrt(frames/3));  % heuristic default

% - may need adjustment!
% Extended-ICA option:
DEFAULT_EXTENDED     = 0;         % default off
DEFAULT_EXTBLOCKS    = 1;         % number of blocks per kurtosis calculation
DEFAULT_NSUB         = 1;         % initial default number of assumed sub-Gaussians for extended-ICA
DEFAULT_EXTMOMENTUM  = 0.5;       % momentum term for computing extended-ICA kurtosis
MAX_KURTSIZE         = 6000;      % max points to use in kurtosis calculation
MIN_KURTSIZE         = 2000;      % minimum good kurtosis size (flag warning)
SIGNCOUNT_THRESHOLD  = 25;        % raise extblocks when sign vector unchanged after this many steps
SIGNCOUNT_STEP       = 2;         % extblocks increment factor

DEFAULT_SPHEREFLAG   = 'on';      % use the sphere matrix as the default
%   starting weight matrix
DEFAULT_PCAFLAG      = 'off';     % don't use PCA reduction
DEFAULT_POSACTFLAG   = 'on';      % use posact()
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule

%--constrained ICA parameters
CONSTRAINED_COMPONENTS=3; % NUMBER OF COMPONENTS FROM EACH DATASET BEING CONSTRAINED
CONSTRAINED_CONNECTION=0.2; % CORRELATION THRESHOLD TO BE CONSTRAINED; HIGH THRESHOLD WILL BE STRENGTHENED.
% ENDURANCE = -1e-2; % the maximumlly allowed descending trend of entropy;
% (Not USing ENDURANCE at all)

%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
epochs = 1;							 % do not care how many epochs in data

pcaflag    = DEFAULT_PCAFLAG;
sphering   = DEFAULT_SPHEREFLAG;     % default flags
posactflag = DEFAULT_POSACTFLAG;
verbose    = DEFAULT_VERBOSE;
block      = DEFAULT_BLOCK;          % heuristic default - may need adjustment!
lrate      = DEFAULT_LRATE;
annealdeg  = DEFAULT_ANNEALDEG;
annealstep = DEFAULT_ANNEALSTEP;     % defaults declared below
nochange   = DEFAULT_STOP;
momentum   = DEFAULT_MOMENTUM;
maxsteps   = DEFAULT_MAXSTEPS;

weights    = 0;                      % defaults defined below
ncomps     = chans;
biasflag   = DEFAULT_BIASFLAG;

extended   = DEFAULT_EXTENDED;
extblocks  = DEFAULT_EXTBLOCKS;
kurtsize   = [MAX_KURTSIZE,MAX_KURTSIZE];
signsbias  = 0.02;                   % bias towards super-Gaussian components
extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
nsub       = DEFAULT_NSUB;
wts_blowup = [0,0,0];                      % flag =1 when weights too large
wts_passed = 0;                      % flag weights passed as argument
%
Connect_threshold =CONSTRAINED_CONNECTION; % set a threshold to select columns constrained.
MaxComCon  =       CONSTRAINED_COMPONENTS;
% *** Endurance is not used ***
% trendPara  = ENDURANCE; %depends on the requirement on connection; the more negative,the stronger the contrains ,that may cause overfitting
%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
if (nargin> 1 & rem(nargin,2) == 0)
    fprintf('runica(): Even number of input arguments???')
    return
end
for i = 3:2:nargin % for each Keyword
    Keyword = eval(['p',int2str((i-3)/2 +1)]);
    Value = eval(['v',int2str((i-3)/2 +1)]);
    if ~isstr(Keyword)
        fprintf('runica(): keywords must be strings')
        return
    end
    Keyword = lower(Keyword); % convert upper or mixed case to lower
    
    if strcmp(Keyword,'weights') | strcmp(Keyword,'weight')
        if isstr(Value)
            fprintf(...
                'runica(): weights value must be a weight matrix or sphere')
            return
        else
            weights = Value;
            wts_passed =1;
        end
    elseif strcmp(Keyword,'ncomps')
        if isstr(Value)
            fprintf('runica(): ncomps value must be an integer')
            return
        end
        if ncomps < urchans & ncomps ~= Value
            fprintf('runica(): Use either PCA or ICA dimension reduction');
            return
        end
        ncomps = Value;
        if ~ncomps,
            ncomps = chans;
        end
    elseif strcmp(Keyword,'pca')
        if ncomps < urchans & ncomps ~= Value
            fprintf('runica(): Use either PCA or ICA dimension reduction');
            return
        end
        if isstr(Value)
            fprintf(...
                'runica(): pca value should be the number of principal components to retain')
            return
        end
        pcaflag = 'on';
        ncomps = Value;
        if ncomps >= chans | ncomps < 1,
            fprintf('runica(): pca value must be in range [1,%d]\n',chans-1)
            return
        end
        chans = ncomps;
    elseif strcmp(Keyword,'posact')
        if ~isstr(Value)
            fprintf('runica(): posact value must be on or off')
            return
        else
            Value = lower(Value);
            if ~strcmp(Value,'on') & ~strcmp(Value,'off'),
                fprintf('runica(): posact value must be on or off')
                return
            end
            posactflag = Value;
        end
    elseif strcmp(Keyword,'lrate')
        if isstr(Value)
            fprintf('runica(): lrate value must be a number')
            return
        end
        lrate = Value;
        if lrate>MAX_LRATE | lrate <0,
            fprintf('runica(): lrate value is out of bounds');
            return
        end
        if ~lrate,
            lrate = DEFAULT_LRATE;
        end
    elseif strcmp(Keyword,'block') | strcmp(Keyword,'blocksize')
        if isstr(Value)
            fprintf('runica(): block size value must be a number')
            return
        end
        block = Value;
        if ~block,
            block = DEFAULT_BLOCK;
        end
    elseif strcmp(Keyword,'stop') | strcmp(Keyword,'nochange') ...
            | strcmp(Keyword,'stopping')
        if isstr(Value)
            fprintf('runica(): stop wchange value must be a number')
            return
        end
        nochange = Value;
    elseif strcmp(Keyword,'maxsteps') | strcmp(Keyword,'steps')
        if isstr(Value)
            fprintf('runica(): maxsteps value must be an integer')
            return
        end
        maxsteps = Value;
        if ~maxsteps,
            maxsteps   = DEFAULT_MAXSTEPS;
        end
        if maxsteps < 0
            fprintf('runica(): maxsteps value must be a positive integer')
            return
        end
    elseif strcmp(Keyword,'anneal') | strcmp(Keyword,'annealstep')
        if isstr(Value)
            fprintf('runica(): anneal step constant must be a number (0,1)')
            return
        end
        annealstep = Value;
        if annealstep <=0 | annealstep > 1,
            fprintf('runica(): anneal step value must be (0,1]')
            return
        end
    elseif strcmp(Keyword,'annealcorr') | strcmp(Keyword,'annealstepcorr')
        if isstr(Value)
            fprintf('runica(): correlation anneal step constant must be a number (0,1)')
            return
        end
        annealcorr = Value;
        if annealcorr <=0 | annealcorr > 1,
            fprintf('runica(): anneal step value must be (0,1]')
            return
        end
    elseif strcmp(Keyword,'annealdeg') | strcmp(Keyword,'degrees')
        if isstr(Value)
            fprintf('runica(): annealdeg value must be a number')
            return
        end
        annealdeg = Value;
        if ~annealdeg,
            annealdeg = DEFAULT_ANNEALDEG;
        elseif annealdeg > 180 | annealdeg < 0
            fprintf('runica(): annealdeg value is out of bounds [0,180]')
            return
            
        end
    elseif strcmp(Keyword,'momentum')
        if isstr(Value)
            fprintf('runica(): momentum value must be a number')
            return
        end
        momentum = Value;
        if momentum > 1.0 | momentum < 0
            fprintf('runica(): momentum value is out of bounds [0,1]')
            return
        end
    elseif strcmp(Keyword,'sphering') | strcmp(Keyword,'sphereing') ...
            | strcmp(Keyword,'sphere')
        if ~isstr(Value)
            fprintf('runica(): sphering value must be on, off, or none')
            return
        else
            Value = lower(Value);
            if ~strcmp(Value,'on') & ~strcmp(Value,'off') & ~strcmp(Value,'none'),
                fprintf('runica(): sphering value must be on or off')
                return
            end
            sphering = Value;
        end
    elseif strcmp(Keyword,'bias')
        if ~isstr(Value)
            fprintf('runica(): bias value must be on or off')
            return
        else
            Value = lower(Value);
            if strcmp(Value,'on')
                biasflag = 1;
            elseif strcmp(Value,'off'),
                biasflag = 0;
            else
                fprintf('runica(): bias value must be on or off')
                return
            end
        end
    elseif strcmp(Keyword,'specgram') | strcmp(Keyword,'spec')
        if isstr(Value)
            fprintf('runica(): specgram argument must be a vector')
            return
        end
        srate = Value(1);
        if (srate < 0)
            fprintf('runica(): specgram srate must be >=0')
            return
        end
        if length(Value)>1
            loHz = Value(2);
            if (loHz < 0 | loHz > srate/2)
                fprintf('runica(): specgram loHz must be >=0 and <= srate/2')
                return
            end
        else
            loHz = 0; % default
        end
        if length(Value)>2
            hiHz = Value(3);
            if (hiHz < loHz | hiHz > srate/2)
                fprintf('runica(): specgram hiHz must be >=loHz and <= srate/2')
                return
            end
        else
            hiHz = srate/2; % default
        end
        if length(Value)>3
            Hzinc = Value(4);
            if (Hzinc<=0 | Hzinc>hiHz-loHz)
                fprintf('runica(): specgram Hzinc must be >0 and <= hiHz-loHz')
                return
            end
        else
            Hzinc = 1; % default
        end
        if length(Value)>4
            Hzframes = Value(5);
            if (Hzframes<0 | Hzframes > size(data1,2))
                fprintf('runica(): specgram frames must be >=0 and <= data length')
                return
            end
        else
            Hzframes = size(data1,2); % default
        end
    elseif strcmp(Keyword,'extended') | strcmp(Keyword,'extend')
        if isstr(Value)
            fprintf('runica(): extended value must be an integer (+/-)')
            return
        else
            extended = 1;      % turn on extended-ICA
            extblocks = fix(Value); % number of blocks per kurt() compute
            if extblocks < 0
                nsub = -1*fix(extblocks);  % fix this many sub-Gauss comps
            elseif ~extblocks,
                extended = 0;             % turn extended-ICA off
            elseif kurtsize>frames,   % length of kurtosis calculation
                kurtsize = frames;
                if kurtsize < MIN_KURTSIZE
                    fprintf(...
                        'runica() warning: kurtosis values inexact for << %d points.\n',...
                        MIN_KURTSIZE);
                end
            end
        end
    elseif strcmp(Keyword,'verbose')
        if ~isstr(Value)
            fprintf('runica(): verbose flag value must be on or off')
            return
        elseif strcmp(Value,'on'),
            verbose = 1;
        elseif strcmp(Value,'off'),
            verbose = 0;
        else
            fprintf('runica(): verbose flag value must be on or off')
            return
        end
    elseif strcmp(Keyword,'dewhitem')
        dewhiteM = Value;
    elseif strcmp(Keyword,'prefs')
        prefs = Value;
    elseif strcmp(Keyword,'tc')
        TC = Value;
    elseif strcmp(Keyword,'whitem')
        whiteM = Value;
    elseif strcmp(Keyword,'constrained_components')
        if(Value>min([size(data1,1),size(data2,1),size(data3,1)] ))
            error(['Data size does not allow more than ',num2str(min([size(data1,1),size(data2,1),size(data3,1)] )),' constrained_components.']);
        end
        MaxComCon  = Value ;
    elseif strcmp(Keyword,'constrained_connection')
        
        Connect_threshold =Value; % set a threshold to select columns constrained.
        
        % *** Endurance is not used ***
        %     elseif strcmp(Keyword,'endurance')
        %         trendPara =Value; %
        
        
    else
        fprintf('runica(): unknown flag')
        return
    end
end
%

if ~exist('dewhiteM','var') | ~exist('whiteM','var')
    help ica_fuse_runica_parallelica2
    return
end

%%%%%%%%%%%%%%%%%%%%%%%% Initialize weights, etc. %%%%%%%%%%%%%%%%%%%%%%%%
%
if ~annealstep,
    if ~extended,
        annealstep = DEFAULT_ANNEALSTEP;     % defaults defined above
    else
        annealstep = DEFAULT_EXTANNEAL;       % defaults defined above
    end
end % else use annealstep from commandline

if ~annealcorr,
    annealcorr = DEFAULT_ANNEALCORR;       % defaults defined above
    
end % else use annealstep from commandline

if ~annealdeg,
    annealdeg  = DEFAULT_ANNEALDEG - momentum*90; % heuristic
    if annealdeg < 0,
        annealdeg = 0;
    end
end

if ncomps >  chans | ncomps < 1
    fprintf('runica(): number of components must be 1 to %d.\n',chans);
    return
end

if weights ~= 0,                    % initialize weights
    % starting weights are being passed to runica() from the commandline
    if verbose,
        fprintf('Using starting weight matrix named in argument list ...\n')
    end
    if  chans>ncomps & weights ~=0,
        [r,c]=size(weights);
        if r~=ncomps | c~=chans,
            fprintf(...
                'runica(): weight matrix must have %d rows, %d columns.\n', ...
                chans,ncomps);
            return;
        end
    end
end;

%%%%%%%%%%%%%%%%%%%%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if frames<chans,
    fprintf('runica(): data length %d < data channels %f!\n',frames,chans)
    return
elseif block < 2,
    fprintf('runica(): block size %d too small!\n',block)
    return
elseif block > frames,
    fprintf('runica(): block size exceeds data length!\n');
    return
elseif floor(epochs) ~= epochs,
    fprintf('runica(): data length is not a multiple of the epoch length!\n');
    return
elseif nsub > ncomps
    fprintf('runica(): there can be at most %d sub-Gaussian components!\n',ncomps);
    return
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf( ...
        '\nModality 1: Input data size [%d,%d] = %d channels, %d frames.', ...
        chans(1),frames(1),chans(1),frames(1));
    fprintf( ...
        '\nModality 2: Input data size [%d,%d] = %d channels, %d frames.', ...
        chans(2),frames(2),chans(2),frames(2));
    fprintf( ...
        '\nModality 3: Input data size [%d,%d] = %d channels, %d frames.\n', ...
        chans(3),frames(3),chans(3),frames(3));
    if strcmp(pcaflag,'on')
        fprintf('After PCA dimension reduction,\n  finding ');
    else
        fprintf('Finding ');
    end
    if ~extended
        fprintf('%d ICA components using logistic ICA.\n',ncomps);
    else % if extended
        fprintf('%d ICA components using extended ICA.\n',ncomps);
        if extblocks > 0
            fprintf(...
                'Kurtosis will be calculated initially every %d blocks using %d data points.\n',...
                extblocks,     kurtsize);
        else
            fprintf(...
                'Kurtosis will not be calculated. Exactly %d sub-Gaussian components assumed.\n',...
                nsub);
        end
    end
    fprintf('Modality 1: Initial learning rate will be %g, block size %d.\n',lrate(1),block(1));
    fprintf('Modality 2: Initial learning rate will be %g, block size %d.\n',lrate(2),block(2));
    fprintf('Modality 3: Initial learning rate will be %g, block size %d.\n',lrate(3),block(3));
    if momentum>0,
        fprintf('Momentum will be %g.\n',momentum);
    end
    fprintf( ...
        'Learning rate will be multiplied by %g whenever angledelta >= %g deg.\n', ...
        annealstep,annealdeg);
    fprintf('Training will end when wchange < %g or after %d steps.\n', ...
        nochange,maxsteps);
    if biasflag,
        fprintf('Online bias adjustment will be used.\n');
    else
        fprintf('Online bias adjustment will not be used.\n');
    end
end

%%%%%%%%%%%%%%%%%%% Perform PCA reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Reducing the data to %d principal dimensions...\n',ncomps);
    [eigenvectors1,eigenvalues1,data1] = pcsquash(data1,ncomps(1)); % changed for three dsatasets
    [eigenvectors2,eigenvalues2,data2] = pcsquash(data2,ncomps(2)); % changed for three datasets
    [eigenvectors3,eigenvalues3,data3] = pcsquash(data3,ncomps(3)); % changed for three datasets
    % make data its projection onto the ncomps-dim principal subspace
end

%%%%%%%%%%%%%%%%%%% Perform sphering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if strcmp(sphering,'on'), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose,
        fprintf('Computing the sphering matrix...\n');
    end
    sphere{1} = 2.0*pinv(sqrtm(cov(data1'))); % find the "sphering" matrix = spher()
    sphere{2} = 2.0*pinv(sqrtm(cov(data2'))); % find the "sphering" matrix = spher()
    sphere{3} = 2.0*pinv(sqrtm(cov(data3')));
    if ~weights,
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
        end
        weights=[];
        weights{1} = eye(ncomps(1),chans(1)); % begin with the identity matrix
        weights{2} = eye(ncomps(2),chans(2)); % begin with the identity matrix
        weights{3} = eye(ncomps(3),chans(3));
    else % weights given on commandline
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
        end
    end
    if verbose,
        fprintf('Sphering the data ...\n');
    end
    
    data1 = sphere{1}*data1;      % actually decorrelate the electrode signals
    data2 = sphere{2}*data2;      % actually decorrelate the electrode signals
    data3 = sphere{3}*data3;
elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~weights
        if verbose,
            fprintf('Using the sphering matrix as the starting weight matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere{1} = 2.0*pinv(sqrtm(cov(data1'))); % find the "sphering" matrix = spher()
        weights{1} = eye(ncomps(1),chans(1))*sphere{1}; % begin with the identity matrix
        sphere{1} = eye(chans(1));                 % return the identity matrix
        sphere{2} = 2.0*pinv(sqrtm(cov(data2'))); % find the "sphering" matrix = spher()
        weights{2} = eye(ncomps(2),chans(2))*sphere{2}; % begin with the identity matrix
        sphere{2} = eye(chans(2));                 % return the identity matrix
        sphere{3} = 2.0*pinv(sqrtm(cov(data3')));
        weights{1} = eye(ncomps(3),chans(3))*sphere{3};
        sphere{3} = eye(chans(2));
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere{1} = eye(chans(1));                 % return the identity matrix
        sphere{2} = eye(chans(2));  % return the identity matrix
        sphere{3} = eye(chans(3));
    end
elseif strcmp(sphering,'none')
    sphere{1} = eye(chans(1));                     % return the identity matrix
    sphere{2} = eye(chans(2));                     % return the identity matrix
    sphere{3} = eye(chans(3));
    if ~weights
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        weights{1} = eye(ncomps(1),chans(1)); % begin with the identity matrix
        weights{2} = eye(ncomps(2),chans(2)); % begin with the identity matrix
        weights{3} = eye(ncomps(3),chans(3));
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
    end
    sphere{1} = eye(chans(1),chans(1));
    sphere{2} = eye(chans(2),chans(2));
    sphere{3} = eye(chans(3),chans(3));
    if verbose,
        fprintf('Returned variable "sphere" will be the identity matrix.\n');
    end
end

%% %%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%  %% Experimental, but wasn't used in the final
% Optimization initialization with three single modalities ICA
% 	disp('---------------------------------------')
% 	disp('Computing initial conditions Modality 1')
% 	disp('---------------------------------------')
% 	[weights{1},sp,~,bias{1},signs{1},t] = icatb_runica(data1, 'verbose','off','lrate',lrate(1));
% 	lrate(1)=t(end);
% 	fprintf('learning rate: %.2e\n',lrate(1))
% 	weights{1} = weights{1}*sp;
% 	disp('---------------------------------------')
% 	disp('Computing initial conditions Modality 2')
%     disp('---------------------------------------')
% 	[weights{2},sp,~,bias{2},signs{2},t] = icatb_runica(data2, 'verbose','off','lrate',lrate(2));
% 	weights{2} = weights{2}*sp;
% 	lrate(2)=t(end);
% 	fprintf('learning rate: %.2e\n',lrate(2))
% 	disp('---------------------------------------')
% 	disp('Computing initial conditions Modality 3')
%     disp('---------------------------------------')
% 	[weights{3},sp,~,bias{3},signs{3},t] = icatb_runica(data3, 'verbose','off','lrate',lrate(3));
% 	weights{3} = weights{3}*sp;
% 	lrate(3)=t(end);
% 	fprintf('learning rate: %.2e\n',lrate(3))

lastt=fix((datalength./block-1).*block+1);
BI{1}=block(1)*eye(ncomps(1),ncomps(1));
BI{2}=block(2)*eye(ncomps(2),ncomps(2));
BI{3}=block(3)*eye(ncomps(3),ncomps(3));
delta{1}=zeros(1,chans(1)*ncomps(1));
delta{2}=zeros(1,chans(2)*ncomps(2));
delta{3}=zeros(1,chans(3)*ncomps(3));
degconst = 180./pi;
startweights = weights;
prevweights = startweights;
oldweights = startweights;
prevwtchange{1} = zeros(chans(1),ncomps(1));
prevwtchange{2} = zeros(chans(2),ncomps(2));
prevwtchange{3} = zeros(chans(3),ncomps(3));
oldwtchange{1} = zeros(chans(1),ncomps(1));
oldwtchange{2} = zeros(chans(2),ncomps(2));
oldwtchange{3} = zeros(chans(3),ncomps(3));

lrates = zeros(2,maxsteps);
onesrow{1} = ones(1,block(1));
onesrow{2} = ones(1,block(2));
onesrow{3} = ones(1,block(3));
bias{1} = zeros(ncomps(1),1);
bias{2} = zeros(ncomps(2),1);
bias{3} = zeros(ncomps(3),1);
signs{1} = ones(1,ncomps(1));    % initialize signs to nsub -1, rest +1
signs{2} = ones(1,ncomps(2));    % initialize signs to nsub -1, rest +1
signs{3} = ones(1,ncomps(3));
for k=1:nsub
    signs{1}(k) = -1;
    signs{2}(k) = -1;
    signs{3}(k) = -1;
end
if extended & extblocks < 0 & verbose,
    fprintf('Fixed extended-ICA sign assignments:  ');
    for k=1:ncomps
        fprintf('%d ',signs(k));
    end; fprintf('\n');
end
signs{1} = diag(signs{1}); % make a diagonal matrix
signs{2} = diag(signs{2}); % make a diagonal matrix
signs{3} = diag(signs{3});
oldsigns{1} = zeros(size(signs{1}));
oldsigns{2} = zeros(size(signs{2}));
oldsigns{3} = zeros(size(signs{3}));

change=[0,0,0];
signcount =[ 0,0,0];              % counter for same-signs
signcounts = [];
urextblocks = extblocks;    % original value, for resets
old_kk{1} = zeros(1,ncomps(1));   % for kurtosis momemtum
old_kk{2} = zeros(1,ncomps(2));   % for kurtosis momemtum
old_kk{3} = zeros(1,ncomps(3));

%%%%%%%% ICA training loop using the logistic sigmoid %%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf('Beginning ICA training ...');
    if extended,
        fprintf(' first training step may be slow ...\n');
    else
        fprintf('\n');
    end
end
step=[0,0,0];
blockno =[1,1,1];  % running block counter for kurtosis interrupts
stopsign=[0,0,0];
angledelta=[0,0,0];
alphak=ones(1,MaxComCon);
Crate=[1,1,1];
ii=1;
sim=[];
iters=0;
num=[0 0 0];
entropychange=[0 0 0];
echange=zeros(540,3);
learnr=zeros(540,3);


% Debug variables
debug_entropy=cell(3,1);
debug_entropy{1} = NaN*zeros(540,2);
debug_entropy{2} = NaN*zeros(540,2);
debug_entropy{3} = NaN*zeros(540,2);

%%
while step(1) < maxsteps || step(2) < maxsteps || step(3) < maxsteps
    iters=iters+1;
    %% -----------------------------------------------------------
    % Find the Loading Matrices
    A1 = WtoA(weights{1},dewhiteM{1});
    A2 = WtoA(weights{2},dewhiteM{2});
    A3 = WtoA(weights{3},dewhiteM{3});
    %% -----------------------------------------------------------
    % Sort the triplets in correlation order
    [idx,simt]=find_mean_v2(A1,A2,A3,MaxComCon); %give me the MaxComCon most corr
    idxt{ii}=idx;
    sim {ii}=simt;
    %% -----------------------------------------------------------
    % Calculate variable for the most correlated triplets
    for trip = 1:MaxComCon
        [sim_c_out,dA_corr1_out,dA_corr2_out,dA_corr3_out,alphak_out] = pre_calc_links_deltas(A1,A2,A3,Crate,idx,trip,Connect_threshold,alphak(trip));
        alphak(trip) = alphak_out;	% alphak is realted to the step length on the optimizatino search
        sim_c(ii,trip,:)=sim_c_out;
        % The Actual Loading Differences
        dA_corr1(:,trip)=dA_corr1_out';
        dA_corr2(:,trip)=dA_corr2_out';
        dA_corr3(:,trip)=dA_corr3_out';
    end
    if verbose
        for trip = 1:MaxComCon
            fprintf('Triplet %d: Correlation between modalities 12 13 and 23: %.2f %.2f %.2f  \n',trip,sim_c(ii,trip,1),sim_c(ii,trip,2),sim_c(ii,trip,3))
        end
    end
    %%
    %  data1 ICA training
    if ~stopsign(1) % && mod(iters,3)==0
        num(1)=num(1)+1;
        
        % the beginning entropy of each step
        entropy(1) = compute_entropy(weights{1},data1, bias{1});
        debug_entropy{1}(num(1),1)=entropy(1);
        
        % begin to update weight matrix
        permuteVec = randperm(datalength(1)); % shuffle data order at each step
        for t=1:block(1):lastt(1), %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
            if biasflag
                
                u=weights{1}*data1(:, permuteVec(t:t+block(1)-1)) + bias{1}*onesrow{1}; % test by jy
                
            else
                u=weights{1}*data1(:, permuteVec(t:t+block(1)-1));
            end
            if ~extended
                %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
                
                y=1./(1+exp(-u));
                weights{1}=weights{1}+lrate(1)*(BI{1}+(1-2*y)*u')*weights{1};                   %
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % extended-ICA
                %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
                y=tanh(u);                                                       %
                weights{1} = weights{1} + lrate(1)*(BI{1}-signs{1}*y*u'-u*u')*weights{1};          %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            if biasflag
                if ~extended
                    %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                    
                    bias{1} = bias{1} + lrate(1)*sum((1-2*y)')'; % for logistic nonlin. %
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else % extended
                    %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    bias{1} = bias{1} + lrate(1)*sum((-2*y)')';  % for tanh() nonlin.   %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            
            if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                weights{1} = weights{1} + momentum*prevwtchange{1};
                prevwtchange{1} = weights{1}-prevweights{1};
                prevweights{1} = weights{1};
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if max(max(abs(weights{1}))) > MAX_WEIGHT
                wts_blowup(1) = 1;
                change(1) = nochange;
            end
            if extended & ~wts_blowup(1)
                %
                %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
                %
                if extblocks > 0 & rem(blockno(1),extblocks) == 0,
                    % recompute signs vector using kurtosis
                    if kurtsize(1) < frames(1)
                        %rp = fix(rand(kurtsize)*datalength);    % pick random subset
                        rp = randperm(datalength(1));
                        % of data. Compute
                        partact=weights{1}*data1(:,rp(1:kurtsize(1))); % partial activation
                    else                                        % for small data sets,
                        partact=weights{1}*data1;                   % use whole data
                    end
                    m2=mean(partact'.^2).^2;
                    m4= mean(partact'.^4);
                    kk= (m4./m2)-3.0;                           % kurtosis estimates
                    if extmomentum
                        kk = extmomentum*old_kk{1} + (1.0-extmomentum)*kk; % use momentum
                        old_kk{1} = kk;
                    end
                    signs{1}=diag(sign(kk+signsbias));             % pick component signs
                    if signs{1} == oldsigns{1},
                        signcount = signcount+1;
                    else
                        signcount = 0;
                    end
                    oldsigns{1} = signs{1};
                    signcounts = [signcounts signcount];
                    if signcount >= SIGNCOUNT_THRESHOLD,
                        extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                        signcount = 0;                             % less frequent if sign
                    end                                         % is not changing
                end % extblocks > 0 & . . .
            end % if extended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blockno(1) = blockno(1) + 1;
            if wts_blowup(1)
                break
            end
        end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---------------------------
        
        % if weight is not  blowup, update
        if ~wts_blowup(1)
            
            oldwtchange{1} = weights{1}-oldweights{1};
            step(1)=step(1)+1;
            lrates(1,step(1)) = lrate(1);
            angledelta(1)=[0];
            delta{1}=reshape(oldwtchange{1},1,chans(1)*ncomps(1));
            change(1)=delta{1}*delta{1}';
            
        end
        %DATA1 blow up restart-------------------------------
        if wts_blowup(1) | isnan(change(1))|isinf(change(1)),  % if weights blow up,
            wts_blowup(1) = 1;
            change(1) = nochange;
            disp('Weights blow up! Modality 1')
            fprintf('');
            step(1) = 0; stopsign(1)=0;               % start again
            change(1) = nochange;
            wts_blowup(1) = 0;                    % re-initialize variables
            blockno(1) = [1];
            lrate(1) = lrate(1)*DEFAULT_RESTART_FAC; % with lower learning rate
            weights{1} = startweights{1};            % and original weight matrix
            oldweights{1} = startweights{1};
            oldwtchange{1} = zeros(chans(1),ncomps(1));
            delta{1}=zeros(1,chans(1)*ncomps(1));
            olddelta = delta;
            extblocks = urextblocks;
            prevweights{1} = startweights{1};
            prevwtchange{1} = zeros(chans(1),ncomps(1));
            lrates(1,:) = zeros(1,maxsteps);
            bias{1} = zeros(ncomps(1),1);
            if extended
                signs{1} = ones(1,ncomps(1));    % initialize signs to nsub -1, rest +1
                for k=1:nsub
                    signs{1}(k) = -1;
                end
                signs{1} = diag(signs{1}); % make a diagonal matrix
                oldsigns{1} = zeros(size(signs{1}));
            end
            if lrate(1)> MIN_LRATE
                r = rank(data1);
                if r<ncomps(1)
                    fprintf('Data has rank %d. Cannot compute %d components.\n',...
                        r,ncomps);
                    return
                else
                    fprintf(...
                        'Lowering learning rate to %g and starting again.\n',lrate);
                end
            else
                fprintf( ...
                    'runica(): QUITTING - weight matrix may not be invertible!\n');
                return;
            end
        else % if DATA1 weights in bounds
            
            
            %testing the trend of entropy term, avoiding the overfitting of correlation
            lossf1(step(1))=compute_entropy(weights{1},data1, bias{1});
            debug_entropy{1}(num(1),2)=lossf1(step(1));
            
            %changes of entropy term added by jingyu
            if step(1)>1
                entropychange(1)=lossf1(step(1))-entropy(1);
            else
                entropychange(1)=1;
            end
            
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step(1)> 2 & ~stopsign(1)
                angledelta(1)=acos((delta{1}*olddelta{1}')/sqrt(change(1)*oldchange(1)));
            end
            if verbose,
                if step(1) > 2,
                    if ~extended,
                        fprintf(...
                            '1:step %d - lrate %.2e, crate %.2e, wchange %.2e, angledelta %4.1f deg, entropy change %.2e\n', ...
                            step(1),lrate(1), Crate(1),change(1),degconst*angledelta(1),entropychange(1));
                    else
                        fprintf(...
                            '1:step %d - lrate %.2e, crate %.2e, wchange %.2e, angledelta %4.1f deg, %d subgauss\n',...
                            step(1),lrate(1), Crate(1),change(1),degconst*angledelta(1),(ncomps(1)-sum(diag(signs{1})))/2);
                    end
                elseif ~extended
                    fprintf(...
                        '1:step %d - lrate %.2e, crate %.2e, wchange %.2e\n',step(1),lrate(1), Crate(1),change(1));
                else
                    fprintf(...
                        '1:step %d - lrate %.2e, crate %.2e, wchange %.2e, %d subgauss\n',...
                        step(1),lrate(1), Crate(1),change(1),(ncomps(1)-sum(diag(signs{1})))/2);
                end % step > 2
            end; % if verbose
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if entropychange(1)<0 %| degconst*angledelta(1) > annealdeg,
                lrate(1) = lrate(1)*annealstep;          % anneal learning rate
                %                 Crate(1) = Crate(1)*annealcorr;
                olddelta{1}   = delta{1};                % accumulate angledelta until
                oldchange(1)  = change(1);               %  annealdeg is reached
            elseif step(1) == 1                     % on first step only
                olddelta{1}   = delta{1};                % initialize
                oldchange(1)  = change(1);
            end
            %
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if step(1) >2 & change(1) < nochange,      % apply stopping rule
                stopsign(1)=1;   % stop when weights stabilize
            elseif step(1)>= maxsteps;
                stopsign(1)=1; % max step
            elseif change(1) > DEFAULT_BLOWUP,      % if weights blow up,
                lrate(1)=lrate(1)*DEFAULT_BLOWUP_FAC;    % keep trying
            end;                                 % with a smaller learning rate
            
            %%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            oldweights{1} = weights{1};
            %added by david
            learnr(num(1),1)=change(1);
            
        end; % end if weights in bounds
        echange(num(1),1)=entropychange(1);
        
        e1 = compute_entropy(weights{1},data1, bias{1});
        %% ************* update based on correlation *****************
        A1 = WtoA(weights{1},dewhiteM{1});
        A1(:,idx(:,1)) = A1(:,idx(:,1)) + dA_corr1;
        weights{1} = AtoW(A1,whiteM{1});
        e2 = compute_entropy(weights{1},data1, bias{1});
        if e2<=e1
            Crate(1)= Crate(1)*annealcorr;
        end
        if e2<=entropy(1)
            lrate(1) =lrate(1)*annealstep;
        end
        
    end
    %%
    %----------------
    %  data2 ICA training
    if ~stopsign(2) %&& mod(iters,3)==1
        num(2)=num(2)+1;
        % the beginning entropy of each step
        entropy(2) = compute_entropy(weights{2},data2, bias{2});
        %------------
        debug_entropy{2}(num(2),1)=entropy(2);
        % begin to update W matrix
        permuteVec = randperm(datalength(2)); % shuffle data order at each step
        for t=1:block(2):lastt(2), %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
            if biasflag
                
                u=weights{2}*data2(:, permuteVec(t:t+block(2)-1)) + bias{2}*onesrow{2};
                
            else
                u=weights{2}*data2(:, permuteVec(t:t+block(2)-1));
            end
            if ~extended
                %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
                
                y=1./(1+exp(-u));
                weights{2}=weights{2}+lrate(2)*(BI{2}+(1-2*y)*u')*weights{2};                   %
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % extended-ICA
                %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
                y=tanh(u);                                                       %
                weights{2} = weights{2} + lrate(2)*(BI{2}-signs{2}*y*u'-u*u')*weights{2};          %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            if biasflag
                if ~extended
                    %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                    
                    bias{2} = bias{2} + lrate(2)*sum((1-2*y)')'; % for logistic nonlin. %
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else % extended
                    %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    bias{2} = bias{2} + lrate(2)*sum((-2*y)')';  % for tanh() nonlin.   %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            
            if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                weights{2} = weights{2} + momentum*prevwtchange{2};
                prevwtchange{2} = weights{2}-prevweights{2};
                prevweights{2} = weights{2};
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if max(max(abs(weights{2}))) > MAX_WEIGHT
                disp('Weights blow up! Modality 2')
                wts_blowup(2) = 1;
                change(2) = nochange;
            end
            if extended & ~wts_blowup(2)
                %
                %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
                %
                if extblocks > 0 & rem(blockno(2),extblocks) == 0,
                    % recompute signs vector using kurtosis
                    if kurtsize(2) < frames(2)
                        %rp = fix(rand(kurtsize)*datalength);    % pick random subset
                        rp = randperm(datalength(2));
                        % of data. Compute
                        partact=weights{2}*data2(:,rp(1:kurtsize(2))); % partial activation
                    else                                        % for small data sets,
                        partact=weights{2}*data2;                   % use whole data
                    end
                    m2=mean(partact'.^2).^2;
                    m4= mean(partact'.^4);
                    kk= (m4./m2)-3.0;                           % kurtosis estimates
                    if extmomentum
                        kk = extmomentum*old_kk{2} + (1.0-extmomentum)*kk; % use momentum
                        old_kk{2} = kk;
                    end
                    signs{2}=diag(sign(kk+signsbias));             % pick component signs
                    if signs{2} == oldsigns{2},
                        signcount = signcount+1;
                    else
                        signcount = 0;
                    end
                    oldsigns{2} = signs{2};
                    signcounts = [signcounts signcount];
                    if signcount >= SIGNCOUNT_THRESHOLD,
                        extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                        signcount = 0;                             % less frequent if sign
                    end                                         % is not changing
                end % extblocks > 0 & . . .
            end % if extended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blockno(2) = blockno(2) + 1;
            if wts_blowup(2)
                break
            end
        end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        debug_entropy{2}(num(2),2)=compute_entropy(weights{2},data2, bias{2});
        
        if max(max(abs(weights{2}))) > MAX_WEIGHT
            wts_blowup(2) = 1;
            change(2) = nochange;
            disp('Weights blow up! Modality 2')
        end
        
        % if weight is not  blowup, update
        if ~wts_blowup(2)
            oldwtchange{2} = weights{2}-oldweights{2};
            step(2)=step(2)+1;
            lrates(2,step(2)) = lrate(2);
            angledelta(2)=[0];
            delta{2}=reshape(oldwtchange{2},1,chans(2)*ncomps(2));
            change(2)=delta{2}*delta{2}'          ;
        end
        
        %DATA2 blow up restart
        if wts_blowup(2) | isnan(change(2)) | isinf(change(2)),  % if weights blow up,
            fprintf('');
            step(2) = 0;     stopsign(2)=0;                % start again
            change(2) = nochange;
            wts_blowup(2) = 0;                    % re-initialize variables
            blockno(2) = [1];
            lrate(2) = lrate(2)*DEFAULT_RESTART_FAC; % with lower learning rate
            weights{2} = startweights{2};            % and original weight matrix
            oldweights{2} = startweights{2};
            oldwtchange{2} = zeros(chans(2),ncomps(2));
            delta{2}=zeros(1,chans(2)*ncomps(2));
            olddelta = delta;
            extblocks = urextblocks;
            prevweights{2} = startweights{2};
            prevwtchange{2} = zeros(chans(2),ncomps(2));
            lrates(2,:) = zeros(1,maxsteps);
            bias{2} = zeros(ncomps(2),1);
            if extended
                signs{2} = ones(1,ncomps(2));    % initialize signs to nsub -1, rest +1
                for k=1:nsub
                    signs{2}(k) = -1;
                end
                signs{2} = diag(signs{2}); % make a diagonal matrix
                oldsigns{2} = zeros(size(signs{2}));
            end
            if lrate(2)> MIN_LRATE
                r = rank(data2);
                if r<ncomps(2)
                    fprintf('Data has rank %d. Cannot compute %d components.\n',...
                        r,ncomps);
                    return
                else
                    fprintf(...
                        'Lowering learning rate to %g and starting again.\n',lrate);
                end
            else
                fprintf( ...
                    'runica(): QUITTING - weight matrix may not be invertible!\n');
                return;
            end
        else % if DATA2 weights in bounds
            
            %testing the trend of entropy term, avoiding the overfitting of correlation
            lossf2(step(2))=compute_entropy(weights{2},data2, bias{2});
            %changes of entropy term added by jingyu
            if step(2)>1;
                entropychange(2)=lossf2(step(2))-entropy(2);
            else
                entropychange(2)=1;
            end
            %--------
            % *** Endurance is not used ***
            %             if entropychange(2)<0
            %                 index1 = ica_fuse_falsemaxdetect(lossf2,trendPara);
            %                 if index1
            %                     Crate(2) = Crate(2)*0.1;         % anneal learning rate empirical
            %                 end
            %             end % end of test------------------------
            
            
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step(2)> 2 & ~stopsign(2)
                angledelta(2)=acos((delta{2}*olddelta{2}')/sqrt(change(2)*oldchange(2)));
            end
            
            
            if verbose,
                if step(2) > 2,
                    if ~extended,
                        fprintf(...
                            '2:step %d - lrate %.2e, crate %.2e, wchange %.2e, angledelta %4.1f deg, entropy change %.2e\n', ...
                            step(2),lrate(2), Crate(2),change(2),degconst*angledelta(2),entropychange(2));
                    else
                        fprintf(...
                            '2:step %d - lrate %.2e, crate %.2e, wchange %.2e, angledelta %4.1f deg, %d subgauss\n',...
                            step(2),lrate(2), Crate(2),change(2),degconst*angledelta(2),(ncomps(2)-sum(diag(signs{2})))/2);
                    end
                elseif ~extended
                    fprintf(...
                        '2:step %d - lrate %.2e, crate %.2e, wchange %.2e\n',step(2),lrate(2), Crate(2),change(2));
                else
                    fprintf(...
                        '2:step %d - lrate %.2e, crate %.2e, wchange %.2e, %d subgauss\n',...
                        step(2),lrate(2), Crate(2),change(2),(ncomps(2)-sum(diag(signs{2})))/2);
                end % step > 2
            end; % if verbose
            %
            
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if entropychange(2)<0 %|degconst*angledelta(2) > annealdeg,
                lrate(2) = lrate(2)*annealstep;          % anneal learning rate
                %                 Crate(2) = Crate(2)*annealcorr;
                olddelta{2}   = delta{2};                % accumulate angledelta until
                oldchange(2)  = change(2);               %  annealdeg is reached
            elseif step(2) == 1                     % on first step only
                olddelta{2}   = delta{2};                % initialize
                oldchange(2)  = change(2);
            end
            %
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if step(2) >2 & change(2) < nochange,      % apply stopping rule
                stopsign(2)=1;               % stop when weights stabilize
            elseif step(2)>= maxsteps
                stopsign(2)=1;               % stop when
            elseif change(2) > DEFAULT_BLOWUP,      % if weights blow up,
                lrate(2)=lrate(2)*DEFAULT_BLOWUP_FAC;    % keep trying
            end;                                 % with a smaller learning rate
            %%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            oldweights{2} = weights{2};
            learnr(num(2),2)=change(2);
            
        end; % end if weights in bounds
        echange(num(2),2)=entropychange(2);
        
        e1 = compute_entropy(weights{2},data2, bias{2});
        %% ********** update based on correlation  ***************
        A2 = WtoA(weights{2},dewhiteM{2});
        A2(:,idx(:,2))=A2(:,idx(:,2)) + dA_corr2;
        weights{2} = AtoW(A2,whiteM{2});
        
        %       A2(:,idx(Cons_com,2))=A2(:,idx(Cons_com,2)) + dA_corr2_1st';
        %       weights{2} = weights{2} + Crate(2)*(AtoW(A2,whiteM{2}) - oldweights{2});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        e2 = compute_entropy(weights{2},data2, bias{2});
        if e2<e1
            Crate(2)= Crate(2)*annealcorr;
        end
        if e2<=entropy(2)
            lrate(2) =lrate(2)*annealstep;
        end
        
        
    end
    %%
    %----------------
    %  data3 ICA training
    if ~stopsign(3) %&& mod(iters,3)==2
        num(3)=num(3)+1;
        % the beginning entropy of each step
        entropy(3)= compute_entropy(weights{3},data3, bias{3});
        debug_entropy{3}(num(3),1)=entropy(3);
        
        %------------
        % begin to update W matrix
        permuteVec = randperm(datalength(3)); % shuffle data order at each step
        for t=1:block(3):lastt(3), %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
            if biasflag
                
                u=weights{3}*data3(:, permuteVec(t:t+block(3)-1)) + bias{3}*onesrow{3};
                
            else
                u=weights{3}*data3(:, permuteVec(t:t+block(3)-1));
            end
            if ~extended
                %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
                
                y=1./(1+exp(-u));
                weights{3}=weights{3}+lrate(3)*(BI{3}+(1-2*y)*u')*weights{3};                   %
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % extended-ICA
                %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
                y=tanh(u);                                                       %
                weights{3} = weights{3} + lrate(3)*(BI{3}-signs{3}*y*u'-u*u')*weights{3};          %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            if biasflag
                if ~extended
                    %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                    
                    bias{3} = bias{3} + lrate(3)*sum((1-2*y)')'; % for logistic nonlin. %
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else % extended
                    %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    bias{3} = bias{3} + lrate(3)*sum((-2*y)')';  % for tanh() nonlin.   %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            
            if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                weights{3} = weights{3} + momentum*prevwtchange{3};
                prevwtchange{3} = weights{3}-prevweights{3};
                prevweights{3} = weights{3};
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if max(max(abs(weights{3}))) > MAX_WEIGHT
                disp('Weights blow up! Modality 3')
                wts_blowup(3) = 1;
                change(3) = nochange;
            end
            if extended & ~wts_blowup(3)
                %
                %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
                %
                if extblocks > 0 & rem(blockno(3),extblocks) == 0,
                    % recompute signs vector using kurtosis
                    if kurtsize(3) < frames(3)
                        %rp = fix(rand(kurtsize)*datalength);    % pick random subset
                        rp = randperm(datalength(3));
                        % of data. Compute
                        partact=weights{3}*data3(:,rp(1:kurtsize(3))); % partial activation
                    else                                        % for small data sets,
                        partact=weights{3}*data3;                   % use whole data
                    end
                    m2=mean(partact'.^2).^2;
                    m4= mean(partact'.^4);
                    kk= (m4./m2)-3.0;                           % kurtosis estimates
                    if extmomentum
                        kk = extmomentum*old_kk{3} + (1.0-extmomentum)*kk; % use momentum
                        old_kk{3} = kk;
                    end
                    signs{3}=diag(sign(kk+signsbias));             % pick component signs
                    if signs{3} == oldsigns{3},
                        signcount = signcount+1;
                    else
                        signcount = 0;
                    end
                    oldsigns{3} = signs{3};
                    signcounts = [signcounts signcount];
                    if signcount >= SIGNCOUNT_THRESHOLD,
                        extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                        signcount = 0;                             % less frequent if sign
                    end                                         % is not changing
                end % extblocks > 0 & . . .
            end % if extended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blockno(3) = blockno(3) + 1;
            if wts_blowup(3)
                break
            end
        end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        debug_entropy{3}(num(3),2)=compute_entropy(weights{3},data3, bias{3});
        
        if max(max(abs(weights{3}))) > MAX_WEIGHT
            wts_blowup(3) = 1;
            change(3) = nochange;
            disp('Weights blow up! Modality 3')
        end
        
        
        % if weight is not  blowup, update
        if ~wts_blowup(3)
            oldwtchange{3} = weights{3}-oldweights{3};
            step(3)=step(3)+1;
            lrates(3,step(3)) = lrate(3);
            angledelta(3)=[0];
            delta{3}=reshape(oldwtchange{3},1,chans(3)*ncomps(3));
            change(3)=delta{3}*delta{3}';%'
        end
        
        %DATA3 blow up restart
        if wts_blowup(3) || isnan(change(3)) || isinf(change(3)),  % if weights blow up,
            fprintf(' Starting Again');
            step(3) = 0;     stopsign(3)=0;                % start again
            change(3) = nochange;
            wts_blowup(3) = 0;                    % re-initialize variables
            blockno(3) = [1];
            lrate(3) = lrate(3)*DEFAULT_RESTART_FAC; % with lower learning rate
            weights{3} = startweights{3};            % and original weight matrix
            oldweights{3} = startweights{3};
            oldwtchange{3} = zeros(chans(3),ncomps(3));
            delta{3}=zeros(1,chans(3)*ncomps(3));
            olddelta = delta;
            extblocks = urextblocks;
            prevweights{3} = startweights{3};
            prevwtchange{3} = zeros(chans(3),ncomps(3));
            lrates(3,:) = zeros(1,maxsteps);
            bias{3} = zeros(ncomps(3),1);
            if extended
                signs{3} = ones(1,ncomps(3));    % initialize signs to nsub -1, rest +1
                for k=1:nsub
                    signs{3}(k) = -1;
                end
                signs{3} = diag(signs{3}); % make a diagonal matrix
                oldsigns{3} = zeros(size(signs{3}));
            end
            if lrate(3)> MIN_LRATE
                r = rank(data3);
                if r<ncomps(3)
                    fprintf('Data has rank %d. Cannot compute %d components.\n',...
                        r,ncomps);
                    return
                else
                    fprintf(...
                        'Lowering learning rate to %g and starting again.\n',lrate);
                end
            else
                fprintf( ...
                    'runica(): QUITTING - weight matrix may not be invertible!\n');
                return;
            end
        else % if DATA3 weights in bounds
            
            %testing the trend of entropy term, avoiding the overfitting of correlation
            lossf3(step(3))=compute_entropy(weights{3},data3, bias{3});
            
            if step(3)>1;
                entropychange(3)=lossf3(step(3))-entropy(3);
            else
                entropychange(3)=1;
            end
            
            
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step(3)> 2 & ~stopsign(3)
                angledelta(3)=acos((delta{3}*olddelta{3}')/sqrt(change(3)*oldchange(3)));
            end
            
            
            if verbose,
                if step(3) > 2,
                    if ~extended,
                        fprintf(...
                            '3:step %d - lrate %.2e, crate %.2e, wchange %.2e, angledelta %4.1f deg, entropy change %.2e \n', ...
                            step(3),lrate(3), Crate(3),change(3),degconst*angledelta(3),entropychange(3));
                    else
                        fprintf(...
                            '3:step %d - lrate %.2e, crate %.2e, wchange %.2e, angledelta %4.1f deg, %d subgauss\n',...
                            step(3),lrate(3), Crate(3),change(3),degconst*angledelta(3),(ncomps(3)-sum(diag(signs{3})))/2);
                    end
                elseif ~extended
                    fprintf(...
                        '3:step %d - lrate %.2e,crate %.2e, wchange %.2e\n',step(3),lrate(3), Crate(3),change(3));
                else
                    fprintf(...
                        '3:step %d - lrate %.2e, crate %.2e, wchange %.2e, %d subgauss\n',...
                        step(3),lrate(3), Crate(3),change(3),(ncomps(3)-sum(diag(signs{3})))/2);
                end % step > 2
            end; % if verbose
            %
            
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if entropychange(3)<0 %|degconst*angledelta(2) > annealdeg,
                lrate(3) = lrate(3)*annealstep;          % anneal learning rate
                %                 Crate(3) = Crate(3)*annealcorr;
                olddelta{3}   = delta{3};                % accumulate angledelta until
                oldchange(3)  = change(3);               %  annealdeg is reached
            elseif step(3) == 1                     % on first step only
                olddelta{3}   = delta{3};                % initialize
                oldchange(3)  = change(3);
            end
            %
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if step(3) >2 & change(3) < nochange,      % apply stopping rule
                stopsign(3)=1;               % stop when weights stabilize
            elseif step(3)>= maxsteps
                stopsign(3)=1;               % stop when
            elseif change(3) > DEFAULT_BLOWUP,      % if weights blow up,
                lrate(3)=lrate(3)*DEFAULT_BLOWUP_FAC;    % keep trying
            end;                                 % with a smaller learning rate
            %%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            oldweights{3} = weights{3};
            learnr(num(3),3)=change(3);
            
        end; % end if weights in bounds
        echange(num(3),3)=entropychange(3);
        
        e1 = compute_entropy(weights{3},data3, bias{3});
        % update based on correlation
        A3 = WtoA(oldweights{3},dewhiteM{3});
        A3(:,idx(:,3))=A3(:,idx(:,3)) + dA_corr3;
        weights{3} = AtoW(A3,whiteM{3});
        
        %           A3(:,idx(Cons_com,3))=A3(:,idx(Cons_com,3)) + dA_corr3_1st';
        %           weights{3} = weights{3} + Crate(3)*(AtoW(A3,whiteM{3}) - oldweights{3});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        e2 = compute_entropy(weights{3},data3, bias{3});
        if e2<=e1
            Crate(3)= Crate(3)*annealcorr;
        end
        if e2<=entropy(3)
            lrate(3) =lrate(3)*annealstep;
        end
        
        
    end
    %%
    %     % Converge if correlations do not change
    %     if ii>1 && abs(mean(sim_c(ii,:))-mean(sim_c(ii-1,:)))<1e-6
    %         disp('Convergence due change on correlation < 1e-6')
    %         stopsign = ones(3,1);
    %     end
    
    %% Plot all Information as it is iterating
    % 		  % Comment if not wanted: plot all correlations among modalities
    % 		  LineW = 5;
    % 		  MarkerS = 15;
    % 		  figure(1);
    % 		  if(iters > 1)
    % 			  hold on;
    % 		  end
    % 		  plot(iters,sim_c(ii,1),'b.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  hold on;
    % 		  grid on;
    % 		  plot(iters,sim_c(ii,2),'r.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  plot(iters,sim_c(ii,3),'m.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  if(iters == 1)
    % 			  title('3way ICA Correlations among Modalities')
    % 			  legend('corr(1,2)','corr(1,3)','corr(2,3)',3);
    % 		  end
    % 		  hold off;
    % 		% Comment if not wanted: plot all entropies
    % 		  figure(2);
    % 		  if(iters > 1)
    % 			  hold on;
    % 		  end
    % 		  plot(iters,compute_entropy(weights{1},data1, bias{1}),'b.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  hold on;
    % 		  grid on;
    % 		  plot(iters,compute_entropy(weights{2},data2, bias{2}),'r.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  plot(iters,compute_entropy(weights{3},data3, bias{3}),'m.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  if(iters == 1)
    % 			  title('3way ICA Entropy for each Modality')
    % 			  legend('H1','H2','H3',3);
    % 		  end
    % 		  hold off;
    % 		% Comment if not wanted: plot all entropies
    % 		  figure(3);
    % 		  if(iters > 1)
    % 			  hold on;
    % 		  end
    % 		  semilogy(iters,change(1),'b.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  hold on;
    % 		  grid on;
    % 		  semilogy(iters,change(2),'r.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  semilogy(iters,change(3),'m.','LineWidth',LineW,'MarkerSize',MarkerS);
    % 		  if(iters == 1)
    % 			  title('3way ICA Weight Change each Modality')
    % 			  legend('WMod1','WMod2','WMod3',3);
    % 		  end
    % 		  hold off;
    % vvOrganizeFigs();
    
    ii=ii+1;
    
    
    if stopsign(1)==1 && stopsign(2)==1 && stopsign(3)==1
        laststep=step;
        step=[maxsteps,maxsteps,maxsteps];                % stop when weights stabilize
    end
    
end %end ica training

if exist('laststep', 'var')
    lrates = lrates(:,1:max(laststep));           % truncate lrate history vector
end
%%%%%%%%%%%%%% Orient components towards positive activation %%%%%%%%%%%
%
if strcmp(posactflag,'on')
    [activations{1},winvout{1},weights{1}] = ica_fuse_posact(data1,weights{1});
    [activations{2},winvout{2},weights{2}] = ica_fuse_posact(data2,weights{2});
    [activations{3},winvout{3},weights{3}] = ica_fuse_posact(data3,weights{3});
    % changes signs of activations and weights to make activations
    % net rms-positive
else
    activations{1} = weights{1}*data1;
    activations{2} = weights{2}*data2;
    activations{3} = weights{3}*data3;
    
end
%%%%%%%%%%%%%% If pcaflag, compose PCA and ICA matrices %%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Composing the eigenvector, weights, and sphere matrices\n');
    fprintf('  into a single rectangular weights matrix; sphere=eye(%d)\n'...
        ,chans);
    weights{1}= weights{1}*sphere{1}*eigenvectors1(:,1:ncomps(1))';
    sphere{1} = eye(urchans(1));
    weights{2}= weights{2}*sphere{2}*eigenvectors2(:,1:ncomps(2))';
    sphere{2} = eye(urchans(2));
    weights{3}= weights{3}*sphere{3}*eigenvectors3(:,1:ncomps(3))';
    sphere{3} = eye(urchans(3));
end
%%%%%% Sort components in descending order of max projected variance %%%%
% DATA1
if verbose,
    fprintf(...
        'Sorting components in descending order of mean projected variance ...\n');
end
if wts_passed == 0
    %
    %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    meanvar  = zeros(ncomps(1),1);      % size of the projections
    if ncomps == urchans % if weights are square . . .
        winv = pinv(weights{1}*sphere{1});
    else
        fprintf('Using pseudo-inverse of weight matrix to rank order component projections.\n');
        winv = pinv(weights{1}*sphere{1});
    end
    for s=1:ncomps(1)
        if verbose,
            fprintf('%d ',s);         % construct single-component data matrix
        end
        % project to scalp, then add row means
        compproj = winv(:,s)*activations{1}(s,:);
        meanvar(s) = mean(sum(compproj.*compproj)/(size(compproj,1)-1));
        % compute mean variance
    end                                         % at all scalp channels
    if verbose,
        fprintf('\n');
    end
    %
    %%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
    %
    [sortvar, windex] = sort(meanvar);
    windex = windex(ncomps(1):-1:1); % order large to small
    meanvar = meanvar(windex);
    %
    %%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
    %
    if nargout>2, % if activations are to be returned
        if verbose,
            fprintf('Permuting the activation wave forms ...\n');
        end
        activations{1} = activations{1}(windex,:);
    else
        clear activations{1}
    end
    weights{1} = weights{1}(windex,:);% reorder the weight matrix
    bias{1}  = bias{1}(windex);		% reorder them
    signs{1} = diag(signs{1});        % vectorize the signs matrix
    signs{1} = signs{1}(windex);      % reorder them
else
    fprintf('Components not ordered by variance.\n');
end
%
%%%%%% Sort components in descending order of max projected variance %%%%
% data 2
if verbose,
    fprintf(...
        'Sorting components in descending order of mean projected variance ...\n');
end
if wts_passed == 0
    %
    %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    meanvar  = zeros(ncomps(2),1);      % size of the projections
    if ncomps == urchans % if weights are square . . .
        winv = pinv(weights{2}*sphere{2});
    else
        fprintf('Using pseudo-inverse of weight matrix to rank order component projections.\n');
        winv = pinv(weights{2}*sphere{2});
    end
    for s=1:ncomps(2)
        if verbose,
            fprintf('%d ',s);         % construct single-component data matrix
        end
        % project to scalp, then add row means
        compproj = winv(:,s)*activations{2}(s,:);
        meanvar(s) = mean(sum(compproj.*compproj)/(size(compproj,1)-1));
        % compute mean variance
    end                                         % at all scalp channels
    if verbose,
        fprintf('\n');
    end
    %
    %%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
    %
    [sortvar, windex] = sort(meanvar);
    windex = windex(ncomps(2):-1:1); % order large to small
    meanvar = meanvar(windex);
    %
    %%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
    %
    if nargout>2, % if activations are to be returned
        if verbose,
            fprintf('Permuting the activation wave forms ...\n');
        end
        activations{2} = activations{2}(windex,:);
    else
        clear activations{2}
    end
    weights{2} = weights{2}(windex,:);% reorder the weight matrix
    bias{2}  = bias{2}(windex);		% reorder them
    signs{2} = diag(signs{2});        % vectorize the signs matrix
    signs{2} = signs{2}(windex);      % reorder them
else
    fprintf('Components not ordered by variance.\n');
end
%%%%%% Sort components in descending order of max projected variance %%%%
% data 3
if verbose,
    fprintf(...
        'Sorting components in descending order of mean projected variance ...\n');
end
if wts_passed == 0
    %
    %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    meanvar  = zeros(ncomps(3),1);      % size of the projections
    if ncomps == urchans % if weights are square . . .
        winv = pinv(weights{3}*sphere{3});
    else
        fprintf('Using pseudo-inverse of weight matrix to rank order component projections.\n');
        winv = pinv(weights{3}*sphere{3});
    end
    for s=1:ncomps(3)
        if verbose,
            fprintf('%d ',s);         % construct single-component data matrix
        end
        % project to scalp, then add row means
        compproj = winv(:,s)*activations{3}(s,:);
        meanvar(s) = mean(sum(compproj.*compproj)/(size(compproj,1)-1));
        % compute mean variance
    end                                         % at all scalp channels
    if verbose,
        fprintf('\n');
    end
    %
    %%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
    %
    [sortvar, windex] = sort(meanvar);
    windex = windex(ncomps(3):-1:1); % order large to small
    meanvar = meanvar(windex);
    %
    %%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
    %
    if nargout>2, % if activations are to be returned
        if verbose,
            fprintf('Permuting the activation wave forms ...\n');
        end
        activations{3} = activations{3}(windex,:);
    else
        clear activations{2}
    end
    weights{3} = weights{3}(windex,:);% reorder the weight matrix
    bias{3}  = bias{3}(windex);		% reorder them
    signs{3} = diag(signs{3});        % vectorize the signs matrix
    signs{3} = signs{3}(windex);      % reorder them
else
    fprintf('Components not ordered by variance.\n');
end
return

end


function alphak=findsteplength1(fk,deltafk,x1,x2,alphak,P,c1,c2)
% 0<c1<0.5<c2<1;
% (f(x+ap)-f(x))/c1/dela_f(x)/P <a
con=1;coml=length(x1);
while con & con<100
    xnew=x1+alphak*P;
    fk1=corrcoef(xnew,x2);fk1=-fk1(2)^2;
    tcov=((xnew-mean(xnew))*(x2-mean(x2))')/(coml-1);
    comterm=2*(tcov)/var(xnew)/var(x2);
    deltafk1=-comterm*(x2-mean(x2)+tcov*(mean(xnew)-xnew)/var(xnew));% 1st order derivative
    
    firstterm1=(fk1-fk)/c1;firstterm2=deltafk*P(:);
    secondterm1=deltafk1*P(:);secondterm2=deltafk*P(:);
    
    if firstterm1 >alphak*firstterm2
        if firstterm1<0;
            alphak = 0.9*abs(firstterm1/firstterm2);
        else
            alphak = 0.1*alphak;
        end
        if alphak<1e-6 alphak=0;con=0;
        else  con=con+1;
        end
        
    elseif secondterm1<0 & secondterm1<c2*secondterm2
        
        % alphak=abs(secondterm1/secondterm2/c2)*alphak;
        con=con+1;
        alphak=1.1*alphak;
    elseif secondterm1 >0 & secondterm2 <0
        alphak=0.9*alphak;    con=con+1;
        
    else
        con=0;
    end
end
end

function [idx,sim]=find_mean_v2(X,Y,Z,num_comps,flag)
% Finds the 'num_comps' most correlated triplet
% Input:
% X,Y,Z: Loadings matrices with size: [num_sub, num_comp]
% num_comps: number of triplets required
% flag (Useless): kept for compatibility with David's function
% Output:
% idx: [num_comps 3] matrix giving the index of the most correlated triplets
% sim: mean correlation of the triplet.
% Author: Alvaro Ulloa
Ac = ica_fuse_combvec(1:size(X,2),1:size(Y,2),1:size(Z,2))'; % all posible triplet combination
nComb = size(Ac,1);% number of all posible triplet combination

c12 = zeros(nComb,1);
c13 = zeros(nComb,1);
c23 = zeros(nComb,1);


rho12 = ica_fuse_corr(X,Y);
rho23 = ica_fuse_corr(Y,Z);
rho13 = ica_fuse_corr(X,Z);

% Listing all possible correlations
for jj = 1:nComb
    c12(jj) = abs(rho12(Ac(jj,1),Ac(jj,2)));
    c13(jj) = abs(rho13(Ac(jj,1),Ac(jj,3)));
    c23(jj) = abs(rho23(Ac(jj,2),Ac(jj,3)));
end

% mean of all possible correlations
sim = abs(mean([c12 c13 c23],2));
[sim_sorted ii_sim_sorted ] =sort(sim,'descend');

idx = Ac(ii_sim_sorted(1:num_comps),:);
sim = sim_sorted(1:num_comps);

end



function entropy = compute_entropy(w, data, bias)
% Computes mean entropy of estimated components. Given the weights (w=A^-1)
% the observed data (X), compute entropy of S=w*X+ bias
% Input:
% w : weight matrix
% data: observed data matrix (subjects by variables)
% bias: bias term for entropy calculation
% Output:
% entropy: mean entropy of components (single value)


u=w*data(:, :) + bias*ones(1,size(data,2));
y=1./(1+exp(-u));
temp=log(abs(w*y.*(1-y))+eps);
entropy=mean(temp(:));
% X = w*data;
% for kk=1:size(X,1)
%     xx = X(kk,:);
%     X(kk,:) = xx./sqrt(xx*xx');
% end
% t = svd(X);
% t = t/sum(t);
% entropy = -sum(t.*log(t));
end

function mx = WtoA(W,dewhiteM)
% Computes the pseudo inverse of W to estimate A in S=AX or WS=X
% Input:
% W: weight matrix
% dewhiteM: un-whitening matrix
% Output:
% mx: Pseudoinverse of W
lambda=0.01;
if cond(W)^-1<2*eps
    %             mx=(dewhiteM{1}*dewhiteM{1}'+sqrt(lambda)*speye(size(dewhiteM{1},1)))';
    [U1,S1,V1]=svd(W);
    S1=diag(lambda*randn(size(W,1),1))+S1;
    W=U1*S1*V1';
    mx = (W' \ dewhiteM')';
else
    mx = (W' \ dewhiteM')';
    %            mx=pinv(dewhiteM{1})'*weights{1};
end

end

function w = AtoW(mx,whiteM)
% inverse function of 'WtoA'

w = whiteM*mx;
w =pinv(w);% weights{1} \ eye(size(weights{1}));

end

% inA1, inA2 and inA3 are the loading matrices
% inCrate
% inidx : list of indexes
% inCons_com : which triplet (1-te most corr, 2 the second most, etc...)
function [sim_d,dA_corr1,dA_corr2,dA_corr3,alphak] = pre_calc_links_deltas(inA1,inA2,inA3,inCrate,inidx,inCons_com,Connect_threshold,inalphak)
sim_d(1)=abs(ica_fuse_corr(inA1(:,inidx(inCons_com,1)),inA2(:,inidx(inCons_com,2))));
sim_d(2)=abs(ica_fuse_corr(inA1(:,inidx(inCons_com,1)),inA3(:,inidx(inCons_com,3))));
sim_d(3)=abs(ica_fuse_corr(inA2(:,inidx(inCons_com,2)),inA3(:,inidx(inCons_com,3))));

c12_flag = sim_d(1)>Connect_threshold;
c13_flag = sim_d(2)>Connect_threshold;
c23_flag = sim_d(3)>Connect_threshold;

a1=inA1(:,inidx(inCons_com,1))';
a2=inA2(:,inidx(inCons_com,2))';
a3=inA3(:,inidx(inCons_com,3))';%

b12=cov(a1,a2); b12=b12(2); tmcorr12=b12/std(a1)/std(a2);
b13=cov(a1,a3); b13=b13(2); tmcorr13=b13/std(a1)/std(a3);
b23=cov(a2,a3); b23=b23(2); tmcorr23=b23/std(a2)/std(a3);

comterm12=2*b12/var(a1)/var(a2);
comterm13=2*b13/var(a1)/var(a3);
comterm23=2*b23/var(a2)/var(a3);

%move 1 to 2
deltaf=comterm12*(a2-mean(a2)+b12*(mean(a1)-a1)/var(a1));% 1st order derivative
P=deltaf./norm(deltaf); %;(H1*deltaf')';
inalphak=findsteplength1(-tmcorr12^2,-deltaf,a1,a2,inalphak,P,0.0001,0.9);
aweights_temp12=inalphak*P; %

%move 1 to 3
deltaf=comterm13*(a3-mean(a3)+b13*(mean(a1)-a1)/var(a1));% 1st order derivative
P=deltaf./norm(deltaf); %;(H1*deltaf')';
inalphak=findsteplength1(-tmcorr13^2,-deltaf,a1,a3,inalphak,P,0.0001,0.9);
aweights_temp13=inalphak*P;

%move 2 to 1
deltaf=comterm12*(a1-mean(a1)+b12*(mean(a2)-a2)/var(a2));% 1st order derivative
P=deltaf./norm(deltaf); %;(H1*deltaf')';
inalphak=findsteplength1(-tmcorr12^2,-deltaf,a2,a1,inalphak,P,0.0001,0.9);
aweights_temp21=inalphak*P;

%move 2 to 3
deltaf=comterm23*(a3-mean(a3)+b23*(mean(a2)-a2)/var(a2));% 1st order derivative
P=deltaf./norm(deltaf); %;(H1*deltaf')';
inalphak=findsteplength1(-tmcorr23^2,-deltaf,a2,a3,inalphak,P,0.0001,0.9);
aweights_temp23=inalphak*P;

%move 3 to 1
deltaf=comterm13*(a1-mean(a1)+b13*(mean(a3)-a3)/var(a3));% 1st order derivative
P=deltaf./norm(deltaf); %;(H1*deltaf')';
inalphak=findsteplength1(-tmcorr13^2,-deltaf,a3,a1,inalphak,P,0.0001,0.9);
aweights_temp31=inalphak*P;%


%move 3 to 2
deltaf=comterm23*(a2-mean(a2)+b23*(mean(a3)-a3)/var(a3));% 1st order derivative
P=deltaf./norm(deltaf); %;(H1*deltaf')';
inalphak=findsteplength1(-tmcorr23^2,-deltaf,a3,a2,inalphak,P,0.0001,0.9);
aweights_temp32=inalphak*P;

dA_corr1=inCrate(1)*(c12_flag*aweights_temp12+ c13_flag*aweights_temp13);
dA_corr2=inCrate(2)*(c12_flag*aweights_temp21+ c23_flag*aweights_temp23);
dA_corr3=inCrate(3)*(c13_flag*aweights_temp31+ c23_flag*aweights_temp32);

alphak=inalphak;
end