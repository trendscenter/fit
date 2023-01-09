% runica() - Perform Independent Component Analysis (ICA) decomposition
%            of psychophysiological data using the infomax ICA algorithm of
%            Bell & Sejnowski (1995) with the natural gradient feature
%            of Amari, Cichocki & Yang, the extended-ICA algorithm
%            of Lee, Girolami & Sejnowski, PCA dimension reduction,
%            and/or specgram() preprocessing (M. Zibulevsky).
% Usage:
%      simply >> [weights,sphere] = runica(data);
%       or
%        else >> [weights,sphere,activations,bias,signs,lrates] ...
%                                 = runica(data,'Key1',Value1',...);
% Input_Variable:
%
% data        = input data (chans,frames*epochs).
%               Note: If data consists of multiple discontinuous epochs,
%               each epoch should be separately baseline-zero'd using:
%                  >> data = rmbase(data,frames,basevector);
%
% Optional_Keywords       Keyword_Values                  Default_Values
%
% 'ncomps'    = [N] number of ICA components to compute (default -> chans)
%               using rectangular ICA decomposition
% 'pca'       = [N] decompose a principal component     (default -> 0=off)
%               subspace of the data. Value is the number of PCs to retain.
% 'sphering'  = ['on'/'off'] flag sphering of data      (default -> 'on')
% 'weights'   = [W] initial weight matrix               (default -> eye())
%                            (Note: if 'sphering' 'off', default -> spher())
% 'lrate'     = [rate] initial ICA learning rate (<< 1) (default -> heuristic)
% 'block'     = [N] ICA block size (<< datalength)      (default -> heuristic)
% 'anneal'    = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended)
%                         controls speed of convergence
% 'annealdeg' = [N] degrees weight change for annealing (default -> 70)
% 'stop'      = [f] stop training when weight-change < this (default -> 1e-6)
% 'maxsteps'  = [N] max number of ICA training steps    (default -> 512)
% 'bias'      = ['on'/'off'] perform bias adjustment    (default -> 'on')
% 'momentum'  = [0<f<1] training momentum               (default -> 0)
% 'extended'  = [N] perform tanh() "extended-ICA" with sign estimation
%               every N training blocks. If N < 0, fix number of sub-Gaussian
%               components to -N [faster than N>0]      (default|0 -> off)
% 'specgram'  = [srate loHz hiHz Hzinc frames] decompose a complex time/frequency
%               transform of the data         (defaults [srate 0 srate/2 1 0])
% 'posact'    = make all component activations net-positive(default 'on'}
% 'verbose'   = give ascii messages ('on'/'off')        (default -> 'on')
%% add by jing sui, 2008.10.21
% 'showscore'= display the weight change, angle change and entropy, default 0
%
% Output_Variables [RO = output in reverse order of projected mean variance
%                        unless starting weight matrix passed ('weights' above)]
%
% weights     = ICA weight matrix (comps,chans)     [RO]
% sphere      = data sphering matrix (chans,chans) = spher(data)
%               Note: unmixing_matrix = weights*sphere {sphering off -> eye(chans)}
% activations = activation time courses of the output components (ncomps,frames*epochs)
% bias        = vector of final (ncomps) online bias [RO]    (default = zeros())
% signs       = extended-ICA signs for components    [RO]    (default = ones())
%                   [-1 = sub-Gaussian; 1 = super-Gaussian]
% lrates      = vector of learning rates used at each training step
%

% Toolbox Citation:
%
% Makeig, Scott et al. "ICA Toolbox for Psychophysiological Research (version 3.4)".
% WWW Site, Computational Neurobiology Laboratory, The Salk Institute for Biological
% Studies <www.cnl.salk.edu/~ica.html>, 1999. [World Wide Web Publication].
%
% 1st Publication:
%
% Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
% "Independent component analysis of electroencephalographic data,"
% In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural
% Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).
%
% For more information:
% http://www.cnl.salk.edu/~scott/icafaq.html - FAQ on ICA/EEG
% http://www.cnl.salk.edu/~scott/icabib.html - mss. on ICA & biosignals
% http://www.cnl.salk.edu/~tony/ica.html - math. mss. on ICA, with kernal code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit history %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  runica()  - by Scott Makeig with contributions from Tony Bell, Te-Won Lee
%              Tzyy-Ping Jung, Sigurd Enghoff, Michael Zibulevsky et al.
%                            CNL / Salk Institute 1996-99
%  04-30-96 built from icatest.m and ~jung/.../wtwpwica.m -sm
%  07-28-97 new runica(), adds bias (default on), momentum (default off),
%           extended-ICA (Lee & Sejnowski, 1997), cumulative angledelta
%           (until lrate drops), keywords, signcount for speeding extended-ICA
%  10-07-97 put acos() outside verbose loop; verbose 'off' wasn't stopping
%  -sm
%  11-11-97 adjusted help msg -sm
%  11-30-97 return eye(chans) if sphering 'off' or 'none' (undocumented option) -sm
%  02-27-98 use pinv() instead of inv() to rank order comps if ncomps < chans -sm
%  04-28-98 added 'posact' and 'pca' flags  -sm
%  07-16-98 reduced length of randperm() for kurtosis subset calc. -se & sm
%  07-19-98 fixed typo in weights def. above -tl & sm
%  12-21-99 added 'specgram' option suggested by Michael Zibulevsky, UNM -sm
%  12-22-99 fixed rand() sizing inefficiency on suggestion of Mike Spratling, UK -sm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

function [weights,sphere,activations,bias,signs,lrates,testrecord, y] = ica_fuse_runica_ccica(data, p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14,p15,v15,p16,v16)
%% make sure the most important vector used first, benefit for the initialization
% data=flipud(data);

%%
if nargin < 1
    help runica
    return
end

[chans frames] = size(data); % determine the data size
urchans = chans;  % remember original data channels
datalength = frames;
%if chans<2 | frames<chans
%   fprintf('\nrunica() - data size (%d,%d) too small.\n\n', chans,frames);
%   return
%end
%
%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      %     anneal by multiplying lrate by this
DEFAULT_EXTANNEAL    = 0.98;      %     or this if extended-ICA
DEFAULT_MAXSTEPS     =512;       % ]top training after this many steps
DEFAULT_MOMENTUM     = 0.0;       % default momentum weight

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.9;       % if weights blowup, restart with lrate
% lower by this factor
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
DEFAULT_LRATE        = 0.015/log(chans);
% heuristic default - may need adjustment
%   for large or tiny data sets!
DEFAULT_BLOCK        = floor(sqrt(frames/3));  % heuristic default
% - may need adjustment!
% Extended-ICA option:
DEFAULT_EXTENDED     = 0;         % default off
DEFAULT_EXTBLOCKS    = 1;         % number of blocks per kurtosis calculation
DEFAULT_NSUB         = 1;         % initial default number of assumed sub-Gaussians
% for extended-ICA
DEFAULT_EXTMOMENTUM  = 0.5;       % momentum term for computing extended-ICA kurtosis
MAX_KURTSIZE         = 6000;      % max points to use in kurtosis calculation
MIN_KURTSIZE         = 2000;      % minimum good kurtosis size (flag warning)
SIGNCOUNT_THRESHOLD  = 25;        % raise extblocks when sign vector unchanged
% after this many steps
SIGNCOUNT_STEP       = 2;         % extblocks increment factor

DEFAULT_SPHEREFLAG   = 'off';      % use the sphere matrix as the default
%   starting weight matrix
DEFAULT_PCAFLAG      = 'off';     % don't use PCA reduction
DEFAULT_POSACTFLAG   = 'on';      % use posact()
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule

DEFAULT_SHOWSCORE=1;           % default, not show weight change and angle change plot.

%% Options for CCICA
% Maximum no. of constraint steps
MAX_CONSTRAINT_STEPS = 15;
% Min correlation between un-mixing matrix before and after adding the
% constraint
MIN_CORRELATION = 0.95;
%the number of components to be contrained
NumConstrainedIC=round(chans/2);
% Note the constrained IC is suggested in range of [NumofIC / 2, NumofIC]
%% End for defining options for CCICA


%
%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargout < 2,
    fprintf('runica() - needs at least two output arguments.\n');
    return
end
epochs = 1;							 % do not care how many epochs in data

pcaflag    = DEFAULT_PCAFLAG;
sphering   = DEFAULT_SPHEREFLAG;     % default flags
posactflag = DEFAULT_POSACTFLAG;
verbose    = DEFAULT_VERBOSE;

block      = DEFAULT_BLOCK;          % heuristic default - may need adjustment!
lrate      = DEFAULT_LRATE;
annealdeg  = DEFAULT_ANNEALDEG;
annealstep = 0;                      % defaults declared below
nochange   = DEFAULT_STOP;
momentum   = DEFAULT_MOMENTUM;
maxsteps   = DEFAULT_MAXSTEPS;

weights    = 0;                      % defaults defined below
ncomps     = chans;
biasflag   = DEFAULT_BIASFLAG;

extended   = DEFAULT_EXTENDED;
extblocks  = DEFAULT_EXTBLOCKS;
kurtsize   = MAX_KURTSIZE;
signsbias  = 0.02;                   % bias towards super-Gaussian components
extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
nsub       = DEFAULT_NSUB;
wts_blowup = 0;                      % flag =1 when weights too large
wts_passed = 1;                      % flag weights passed as argument
testrecord=0;

showscore=DEFAULT_SHOWSCORE;      %% DEFAULT_SHOWSCORE =0,  ADDED BY JING SUI


max_constraint_steps = MAX_CONSTRAINT_STEPS;
min_corr_constraint = MIN_CORRELATION;

%
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
        if length(data)>10000
            nochange=0.000003;
        else
            nochange   = Value;
        end
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
            if (Hzframes<0 | Hzframes > size(data,2))
                fprintf('runica(): specgram frames must be >=0 and <= data length')
                return
            end
        else
            Hzframes = size(data,2); % default
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
    elseif strcmp(Keyword,'showscore')
        if isstr(Value)
            fprintf('runica(): showscore must be 1 or zero')
            return
        else
            showscore=Value;
        end
    elseif strcmpi(Keyword, 'dewhitem')
        %% Dewhitening matrix
        dewhiteM = Value;
    elseif strcmpi(Keyword, 'num_subjects')
        %% Number of subjects
        nsubj = Value;
    elseif strcmpi(Keyword, 'max_constraint_steps')
        max_constraint_steps = Value;
    elseif strcmpi(Keyword, 'min_correlation')
        min_corr_constraint = Value;
    elseif strcmpi(Keyword, 'num_constrained_components')
        NumConstrainedIC = Value;
    else
        fprintf('runica(): unknown flag')
        return
    end
end

%% Error check for dewhitening matrix and groups vector
if ~exist('dewhiteM', 'var')
    error('Dewhitening matrix is not passed as a parameter');
end

if ~exist('nsubj', 'var')
    error('Groups vector containing number of subjects of each group is not passed');
end

%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize weights, etc. %%%%%%%%%%%%%%%%%%%%%%%%
%
if ~annealstep,
    if ~extended,
        annealstep = DEFAULT_ANNEALSTEP;     % defaults defined above
    else
        annealstep = DEFAULT_EXTANNEAL;       % defaults defined above
    end
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
%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf( ...
        '\nInput data size [%d,%d] = %d channels, %d frames.\n', ...
        chans,frames,chans,frames);
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
    fprintf('Initial learning rate will be %g, block size %d.\n',lrate,block);
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
%
%%%%%%%%%%%%%%%%%%%%%%%%% Remove overall row means %%%%%%%%%%%%%%%%%%%%%%%%
%
% if verbose,
%    fprintf('Removing mean of each channel ...\n');
% end
disp('Not removing mean of each channel!!!');
%data = data - mean(data')'*ones(1,frames);      % subtract row means

if verbose,
    fprintf('Final training data range: %g to %g\n', ...
        min(min(data)),max(max(data)));
end
%
%%%%%%%%%%%%%%%%%%% Perform PCA reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Reducing the data to %d principal dimensions...\n',ncomps);
    [eigenvectors,eigenvalues,data] = pcsquash(data,ncomps);
    % make data its projection onto the ncomps-dim principal subspace
end
%
%%%%%%%%%%%%%%%%%%% Perform specgram transformation %%%%%%%%%%%%%%%%%%%%%%%
%
if exist('srate')
    % [P F T] = SPECGRAM(A,NFFT,Fs,WINDOW,NOVERLAP)
    Hzwinlen =  fix(srate/Hzinc);
    Hzfftlen = 2^(ceil(log(Hzwinlen)/log(2)));
    if (Hzwinlen>Hzframes)
        Hzwinlen = Hzframes;
    end
    Hzoverlap = 0;
    if (Hzwinlen>Hzframes)
        Hzoverlap = Hzframes - Hzwinlen;
    end
    % get freqs and times
    [tmp,freqs,tms] = specgram(data(1,:),Hzfftlen,srate,Hzwinlen,Hzoverlap);
    fs = find(freqs>=loHz & freqs <= hiHz);
    % fprintf('   size(fs) = %d,%d\n',size(fs,1),size(fs,2));
    % fprintf('   size(tmp) = %d,%d\n',size(tmp,1),size(tmp,2));
    specdata = reshape(tmp(fs,:),1,length(fs)*size(tmp,2));
    specdata = [real(specdata) imag(specdata)];
    for ch=2:chans
        [tmp] = specgram(data(ch,:),Hzwinlen,srate,Hzwinlen,Hzoverlap);
        tmp = reshape((tmp(fs,:)),1,length(fs)*size(tmp,2));
        specdata = [specdata;[real(tmp) imag(tmp)]]; % channels are rows
    end
    fprintf('Converted data to %d channels by %d=2*%dx%d points spectrogram data.\n',chans,2*length(fs)*length(tms),length(fs),length(tms));
    fprintf('   Low Hz %g, high Hz %g, Hz inc %g, window length %d\n',freqs(fs(1)),freqs(fs(end)),freqs(fs(2))-freqs(fs(1)),Hzwinlen);
    data = specdata;
    datalength=size(data,2);
end
%
%%%%%%%%%%%%%%%%%%% Perform sphering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if strcmp(sphering,'on'), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose,
        fprintf('Computing the sphering matrix...\n');
    end
    sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
    if ~weights,
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
        end
        weights = eye(ncomps,chans); % begin with the identity matrix
    else % weights given on commandline
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
        end
    end
    if verbose,
        fprintf('Sphering the data ...\n');
    end
    data = sphere*data;      % actually decorrelate the electrode signals
    
elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~weights
        if verbose,
            fprintf('Using the sphering matrix as the starting weight matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
        weights = eye(ncomps,chans)*sphere; % begin with the identity matrix
        sphere = eye(chans);                 % return the identity matrix
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere = eye(chans);                 % return the identity matrix
    end
elseif strcmp(sphering,'none')
    sphere = eye(chans);                     % return the identity matrix
    if ~weights
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        weights = eye(ncomps,chans); % begin with the identity matrix
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
    end
    sphere = eye(chans,chans);
    if verbose,
        fprintf('Returned variable "sphere" will be the identity matrix.\n');
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%
lastt=fix((datalength/block-1)*block+1);
BI=block*eye(ncomps,ncomps);
delta=zeros(1,chans*ncomps);
changes = [];
angle=[];
degconst = 180./pi;
startweights = weights;
prevweights = startweights;
oldweights = startweights;
prevwtchange = zeros(chans,ncomps);
oldwtchange = zeros(chans,ncomps);
lrates = zeros(1,maxsteps);
onesrow = ones(1,block);
bias = zeros(ncomps,1);
signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
for k=1:nsub
    signs(k) = -1;
end
if extended & extblocks < 0 & verbose,
    fprintf('Fixed extended-ICA sign assignments:  ');
    for k=1:ncomps
        fprintf('%d ',signs(k));
    end; fprintf('\n');
end
signs = diag(signs); % make a diagonal matrix
oldsigns = zeros(size(signs));;
signcount = 0;              % counter for same-signs
signcounts = [];
urextblocks = extblocks;    % original value, for resets
old_kk = zeros(1,ncomps);   % for kurtosis momemtum
%
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
step=0;
laststep=0;
blockno = 1;  % running block counter for kurtosis interrupts
wmatrix=zeros(chans); % save weights
blowtime=0;
%%
while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    permute=randperm(datalength); % shuffle data order at each step
    
    for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
        if biasflag
            u=weights*data(:,permute(t:t+block-1)) + bias*onesrow;
        else
            u=weights*data(:,permute(t:t+block-1));
        end
        if ~extended
            %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
            y=1./(1+exp(-u));                                                %
            weights=weights+lrate*(BI+(1-2*y)*u')*weights;                   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else % extended-ICA
            %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
            y=tanh(u);                                                       %
            weights = weights + lrate*(BI-signs*y*u'-u*u')*weights;          %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if biasflag
            if ~extended
                %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                bias = bias + lrate*sum((1-2*y)')'; % for logistic nonlin. %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % extended
                %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                bias = bias + lrate*sum((-2*y)')';  % for tanh() nonlin.   %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        
        if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            weights = weights + momentum*prevwtchange;
            prevwtchange = weights-prevweights;
            prevweights = weights;
        end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if max(max(abs(weights))) > MAX_WEIGHT
            wts_blowup = 1;
            change = nochange;
        end
        if extended & ~wts_blowup
            %
            %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
            %
            if extblocks > 0 & rem(blockno,extblocks) == 0,
                % recompute signs vector using kurtosis
                if kurtsize < frames
                    %rp = fix(rand(kurtsize)*datalength);    % pick random subset
                    rp = randperm(datalength);
                    % of data. Compute
                    partact=weights*data(:,rp(1:kurtsize)); % partial activation
                else                                        % for small data sets,
                    partact=weights*data;                   % use whole data
                end
                m2=mean(partact'.^2).^2;
                m4= mean(partact'.^4);
                kk= (m4./m2)-3.0;                           % kurtosis estimates
                if extmomentum
                    kk = extmomentum*old_kk + (1.0-extmomentum)*kk; % use momentum
                    old_kk = kk;
                end
                signs=diag(sign(kk+signsbias));             % pick component signs
                if signs == oldsigns,
                    signcount = signcount+1;
                else
                    signcount = 0;
                end
                oldsigns = signs;
                signcounts = [signcounts signcount];
                if signcount >= SIGNCOUNT_THRESHOLD,
                    extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                    signcount = 0;                             % less frequent if sign
                end                                         % is not changing
            end % extblocks > 0 & . . .
        end % if extended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        blockno = blockno + 1;
        if wts_blowup
            break
        end
    end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   load('IM5_pca_comb_1.mat');
    if ~wts_blowup
        if step<=max_constraint_steps
            u=weights*data;
            y=1./(1+exp(u));
            temp=log(abs(weights*(1-y)*y'));
            enpy(step+1)=mean(temp(:));
        end
        %%       determine lrate at step 8
        if step==max_constraint_steps
            % load('realIM8_pca_comb_1.mat');
            %nochange=0.000001;
            %             load('CC16_pca_comb_1.mat');
            % load('PCAfile.mat');
            %             nsubj=size(whiteM,2);
            %             % determine sh and sp
            %             sh=nsubj/2;
            %             sp=nsubj-sh;
            sh=nsubj(1);
            sp=nsubj(2);
            nsubject=sh+sp;
            whiteM=pinv(dewhiteM);
            %            sh=round(size(chans,1)/2);
            %            sp=size(chans,1)-sh;
            %            nsubject=sh+sp;
            %            dewhiteM=eye(size(chans,1));
            %            whiteM=eye(size(chans,1));
            %% training lrate
            initial_lrate=0.1/log(chans)*log10(length(data));
            lrstep=1;
            min_lr=0.0005;
            MAX_W_trainLR=10000;
            while initial_lrate(lrstep)>min_lr
                %% when resume training learningrate, weights should use the original one
                w=weights;
                blockno = 1;
                if biasflag
                    b=bias;
                end
                lrate2=initial_lrate(lrstep);
                wts_blowup=0;
                for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
                    if biasflag
                        u=w*data(:,permute(t:t+block-1)) + b*onesrow;
                    else
                        u=w*data(:,permute(t:t+block-1));
                    end
                    if ~extended
                        %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%
                        y=1./(1+exp(-u));
                        w=w+lrate*(BI+(1-2*y)*u')*w;
                        AA=dewhiteM*inv(w);
                        A11=AA(1:sh,:);
                        A22=AA(sh+1:end,:);
                        %% sort by p value of ttest, decide which component to be constrained
                        fenmu=ica_fuse_nanvar(A11)*(sh-1)+ica_fuse_nanvar(A22)*(sp-1);
                        fenzi=(mean(A11)-mean(A22)).^2;
                        T1=fenzi./fenmu;
                        [t,ind]=sort(T1);
                        indt=ind(end:-1:end-NumConstrainedIC+1);
                        pbefore(step,:)=2*ica_fuse_spm_Tcdf(-sqrt(T1(indt)*sh*sp*(nsubject-2)/nsubject),nsubject-2);
                        Tsumsqr=sum(T1(indt)*sh*sp*(nsubject-2)/nsubject);
                        %% Update A1A2 by current learning rate
                        for i=indt
                            a1=A11(:,i);
                            a2=A22(:,i);
                            commonterm=(mean(a1)-mean(a2))/(ica_fuse_nanvar(a1)*(sh-1)+ica_fuse_nanvar(a2)*(sp-1));
                            deltafk1=2*(-commonterm/sh+commonterm.^2*(a1-mean(a1)));
                            deltafk2=2*(commonterm/sp+commonterm.^2*(a2-mean(a2)));
                            a1=a1-lrate2*deltafk1;
                            a2=a2-lrate2*deltafk2;
                            A11(:,i)=a1;
                            A22(:,i)=a2;
                        end
                        AA=[A11;A22];
                        w=inv(whiteM*AA);
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else % extended-ICA
                        %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
                        y=tanh(u);                                                       %
                        w = w + lrate*(BI-signs*y*u'-u*u')*w;          %
                        %% test learning rate
                        AA=dewhiteM*inv(w);
                        A11=AA(1:sh,:);
                        A22=AA(sh+1:end,:);
                        %% sort by p value of ttest, decide which component to be constrained
                        fenmu=ica_fuse_nanvar(A11)*(sh-1)+ica_fuse_nanvar(A22)*(sp-1);
                        fenzi=(mean(A11)-mean(A22)).^2;
                        T1=fenzi./fenmu;
                        [t,ind]=sort(T1);
                        indt=ind(end:-1:end-NumConstrainedIC+1);
                        Tsumsqr=sum(T1(indt)*sh*sp*(nsubject-2)/nsubject);
                        pbefore(step,:)=2*ica_fuse_spm_Tcdf(-sqrt(T1(indt)*sh*sp*(nsubject-2)/nsubject),nsubject-2);
                        %% Update A1A2 by current learning rate
                        for i=indt
                            a1=A11(:,i);
                            a2=A22(:,i);
                            commonterm=(mean(a1)-mean(a2))/(ica_fuse_nanvar(a1)*(sh-1)+ica_fuse_nanvar(a2)*(sp-1));
                            deltafk1=2*(-commonterm/sh+commonterm.^2*(a1-mean(a1)));
                            deltafk2=2*(commonterm/sp+commonterm.^2*(a2-mean(a2)));
                            a1=a1-lrate2*deltafk1;
                            a2=a2-lrate2*deltafk2;
                            A11(:,i)=a1;
                            A22(:,i)=a2;
                        end
                        AA=[A11;A22];
                        w=inv(whiteM*AA);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    if biasflag
                        if ~extended
                            %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                            b = b + lrate*sum((1-2*y)')'; % for logistic nonlin. %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else % extended
                            %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            b = b + lrate*sum((-2*y)')';  % for tanh() nonlin.   %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
                    end
                    corrw(lrstep)=abs(ica_fuse_corr2(w,weights));
                    
                    if (corrw < min_corr_constraint)
                        wts_blowup = 1;
                        change = nochange;
                    end
                    blockno = blockno + 1;
                    if wts_blowup
                        break
                    end
                end
                if wts_blowup
                    lrate2=initial_lrate(lrstep)*0.9;
                    lrstep=lrstep+1;
                    initial_lrate(lrstep)=lrate2;
                    if lrate2<min_lr
                        blowtime=blowtime+1;
                        close all;
                        if blowtime>=2
                            disp('use minimum learningrate for CCICA');
                            lrate2=min_lr;
                            wts_blowup=0;
                            if showscore
                                scrsz = get(0,'ScreenSize');
                                figure('Position',[scrsz(3)/8 scrsz(4)/16 scrsz(3)*5/8 scrsz(4)*7/8]);
                                set(gcf,'defaultaxesfontsize',12);
                                set(gcf,'defaultaxesfontweight','bold');
                            end
                            
                            break;
                        else
                            % step=step-1;
                            wts_blowup=1;
                            
                            break;
                        end
                    end
                end
                u=w*data;
                y=1./(1+exp(u));
                temp=log(abs(w*(1-y)*y'));
                enpy1(lrstep)=mean(temp(:));
                u=weights*data;
                y=1./(1+exp(u));
                temp=log(abs(weights*(1-y)*y'));
                enpy2=mean(temp(:));
                if ~wts_blowup
                    lrate2=initial_lrate(lrstep);
                    fprintf('\n');
                    disp('Add contraint on mixing coeffcients') ;
                    fprintf('Learningrate= %5.4f; Correlation between Wim and W = %5.4f. \n',lrate2,corrw(lrstep));
                    disp('Sorting components via T2 in descending order, index are: ');
                    disp(indt);
                    fprintf('\n')
                    %                         if enpy1>=enpy2*0.9
                    %                        [enpy1(lrstep),enpy2,lrate2],corrw(lrstep),indt,
                    break
                    %                     else
                    %                         %  [enpy1,enpy2,lrate2]
                    %                         lrate2=initial_lrate(lrstep)*0.9;
                    %                         lrstep=lrstep+1;
                    %                         initial_lrate(lrstep)=lrate2;
                    %                     end
                    %                                         [enpy1,enpy2]
                    %                                         lrate2,corr2(w,weights),
                    %                                         break
                end
                %%
                if showscore
                    if lrstep>1
                        if lrstep==2
                            scrsz = get(0,'ScreenSize');
                            figure('Position',[scrsz(3)/8 scrsz(4)/32 scrsz(3)*5/8 scrsz(4)*7/8]);
                            set(gcf,'defaultaxesfontsize',12);
                            set(gcf,'defaultaxesfontweight','bold');
                            grid on;
                        end
                        axes('position',[0.02, 0.66, 0.95 0.33]);
                        subplot(311),plot(corrw,'r-*'),
                        title('weights correlation');ylabel('abs corr');xlabel('step');axis tight;%text(0,2,num2str(corrw(end)));
                        axes('position',[0.02, 0.33, 0.95 0.33]);
                        subplot(312),plot(initial_lrate,'b--.'),title('constraint strength');xlabel('step');axis tight;
                        axes('position',[0.02, 0, 0.95 0.33]);
                        subplot(313),plot(enpy1,'m:.'),title('entropy from updated w' );xlabel('step');axis tight;
                        hold on,plot(enpy2*ones(1,lrstep),'g-');
                        drawnow
                    end
                end
                %%
            end
            
            %% training end
            %    lrate2=0.16;
            learningrate(1:step)=lrate2;
        end
        %% after step 8, use CCICA
        if step>max_constraint_steps
            %            epwindow=enpy(step-8:step);
            %%
            wchwindow=wchange(step-8:step);
            x=[ones(1,9);1:9]';
            aa=x\wchwindow';
            trendrecord(step-8)=aa(2);
            %%
            w=weights;
            AA=dewhiteM*inv(w);
            A1=AA(1:sh,:);
            A2=AA(sh+1:end,:);
            %% vartype=???
            fenmu=ica_fuse_nanvar(A1)*(sh-1)+ica_fuse_nanvar(A2)*(sp-1);
            fenzi=(mean(A1)-mean(A2)).^2;
            T1=fenzi./fenmu;
            
            % sort by p value of ttest2
            [p,ind]=sort(T1);
            indt=ind(end:-1:end-NumConstrainedIC+1);
            % indt=ind(end);
            pbefore(step,:)=2*ica_fuse_spm_Tcdf(-sqrt(T1(indt)*sh*sp*(nsubject-2)/nsubject),nsubject-2);
            Tsumsqr_bf(step)=sum(T1(indt)*sh*sp*(nsubject-2)/nsubject);
            % test p value of ttest again
            % ind=find(learningrate);
            lrate2=learningrate(end);
            % [lrate2,A1,A2,eptemp]=findrate(lrate2,A1,A2,data,whiteM,epwindow,indt);
            [lrate2,A1,A2,change]=ica_fuse_findratebyw(lrate2,A1,A2,oldweights,wchwindow,whiteM,indt);
            learningrate=[learningrate,lrate2];
            %           enpy(step+1)=eptemp;
            %             add constraint to A
            fenmu=ica_fuse_nanvar(A1)*(sh-1)+ica_fuse_nanvar(A2)*(sp-1);
            fenzi=(mean(A1)-mean(A2)).^2;
            T2=fenzi./fenmu*sh*sp*(nsubject-2)/nsubject;
            
            pafter(step,:)=2*ica_fuse_spm_Tcdf(-sqrt(T2(indt)),nsubject-2);
            %             test divergence between groups
            %            reture to weights
            AA=[A1;A2];
            w=inv(whiteM*AA);
            u=w*data;
            y=1./(1+exp(u));
            temp=log(abs(w*(1-y)*y'));
            enpy(step)=mean(temp(:));
            Tsumsqr(step)=sum(T2(indt));
        end
        
        oldwtchange = weights-oldweights;
        step=step+1;
        %
        %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
        %
        lrates(step) = lrate;
        angledelta=0.;
        delta=reshape(oldwtchange,1,chans*ncomps);
        change=delta*delta';
        wchange(step+1)=change;
        wmatrix(:,:,step)=weights;
    end
    %
    %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
    %
    if wts_blowup | isnan(change)|isinf(change),  % if weights blow up,
        fprintf('');
        step = 0;                          % start again
        change = nochange;
        wts_blowup = 0;                    % re-initialize variables
        blockno = 1;
        lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
        weights = startweights;            % and original weight matrix
        oldweights = startweights;
        change = nochange;
        oldwtchange = zeros(chans,ncomps);
        delta=zeros(1,chans*ncomps);
        olddelta = delta;
        extblocks = urextblocks;
        prevweights = startweights;
        prevwtchange = zeros(chans,ncomps);
        lrates = zeros(1,maxsteps);
        bias = zeros(ncomps,1);
        if extended
            signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
            for k=1:nsub
                signs(k) = -1;
            end
            signs = diag(signs); % make a diagonal matrix
            oldsigns = zeros(size(signs));;
        end
        if lrate> MIN_LRATE
            r = rank(data);
            if r<ncomps
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
    else % if weights in bounds
        %
        %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
        %
        if step> 2
            angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
            angle=[angle,angledelta];
            if degconst*angledelta<60
                if step>max_constraint_steps+1
                    learningrate(end)=lrate2*0.95;
                end
            end
        end
        if verbose,
            if step > 2,
                if ~extended,
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
                        step,lrate,change,degconst*angledelta);
                else
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\n',...
                        step,lrate,change,degconst*angledelta,(ncomps-sum(diag(signs)))/2);
                end
            elseif ~extended
                fprintf(...
                    'step %d - lrate %5f, wchange %7.6f\n',step,lrate,change);
            else
                fprintf(...
                    'step %d - lrate %5f, wchange %7.6f, %d subgauss\n',...
                    step,lrate,change,(ncomps-sum(diag(signs)))/2);
            end % step > 2
        end; % if verbose
        %
        %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        changes = [changes change];
        oldweights = weights;
        %
        %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if degconst*angledelta > annealdeg,
            lrate = lrate*annealstep;          % anneal learning rate
            olddelta   = delta;                % accumulate angledelta until
            oldchange  = change;               %  annealdeg is reached
        elseif step == 1                     % on first step only
            olddelta   = delta;                % initialize
            oldchange  = change;
        end
        %% plot variable trend
        if showscore
            if step>=max_constraint_steps
                %           hold on,
                %           axes('position',[0.02, 0.50, 0.95 0.45]);
                %           subplot(211),plot(step,log10(change),'r-s'),title('weights change');axis tight;
                %           hold on,
                %           axes('position',[0.02, 0.01, 0.95 0.45]);
                %           subplot(212),plot(step,angledelta*degconst,'b--o'),title('delta angle');axis tight;
                %           drawnow
                if step==max_constraint_steps
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[scrsz(3)/8 scrsz(4)/16 scrsz(3)*5/8 scrsz(4)*7/8]);
                    set(gcf,'defaultaxesfontsize',12);
                    set(gcf,'defaultaxesfontweight','bold');
                    grid on;
                end
                axes('position',[0.02, 0.75, 0.95 0.24]);
                subplot(411),plot(log10(changes),'r-*'),
                hold on,plot(log10(lrates(1:step)),'g--.'),
                if step>max_constraint_steps
                    hold on,plot(log10(learningrate),'b-.');%legend('cons strength');
                end
                title('weights change');legend('w change','learning rate');ylabel('log10(value)');xlabel('step');axis tight;
                axes('position',[0.02, 0.5, 0.95 0.24]);
                subplot(412),plot(angle*degconst,'b--.'),title('delta angle');ylabel('degrees');xlabel('step');axis tight;
                axes('position',[0.02, 0.25, 0.95 0.24]);
                subplot(413),plot(enpy,'m:.'),title('entropy of CCICA cons6' );xlabel('step');axis tight;
                drawnow
                if step>16
                    axes('position',[0.02, 0, 0.95 0.24]);
                    subplot(414),plot(Tsumsqr(16:end),'k-.'),
                    hold on, plot(Tsumsqr_bf(16:end),'r:.');
                    title('square T');xlabel('step');axis tight;
                    drawnow;
                end
            end
        end
        %%
        %
        %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if step >2 & change < nochange,      % apply stopping rule
            change=0.000001;
            laststep=step;
            step=maxsteps;                  % stop when weights stabilize
        elseif change > DEFAULT_BLOWUP,      % if weights blow up,
            lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying
        end;                                 % with a smaller learning rate
    end; % end if weights in bounds
    
end;
% end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if ~laststep
    laststep = step;
end;
lrates = lrates(1,1:laststep);           % truncate lrate history vector
%% plot constraints
% figure,plot(enpy);
% title(['entropy of CCICA cons6 at ini_lr=',num2str(lrate2)]);
% ind=find(enpy==max(enpy));
% enpy(ind),ind
% figure,plot(log10(changes),'r');
% title('log changes')
% figure,
% for i=1:length(indt)
%     hold on, plot(pbefore(:,i)+i*2,'b');
%     hold on, plot(pafter(:,i)+i*2,'r');
% end
% figure,plot(learningrate);
% title('learning rate');
save('constraint_record.mat','enpy','changes','wmatrix','pbefore','pafter','learningrate','indt');

% figure,plot(trendrecord);
% title('change trend every ten steps');
%
%%%%%%%%%%%%%% Orient components towards positive activation %%%%%%%%%%%
%
if strcmp(posactflag,'on')
    [activations, winvout, weights] = ica_fuse_posact(data, weights);
    % changes signs of activations and weights to make activations
    % net rms-positive
else
    activations = weights*data;
end
%
%%%%%%%%%%%%%% If pcaflag, compose PCA and ICA matrices %%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Composing the eigenvector, weights, and sphere matrices\n');
    fprintf('  into a single rectangular weights matrix; sphere=eye(%d)\n'...
        ,chans);
    weights= weights*sphere*eigenvectors(:,1:ncomps)';
    sphere = eye(urchans);
end
%
%%%%%% Sort components in descending order of max projected variance %%%%
%
if verbose,
    fprintf(...
        'Sorting components in descending order of mean projected variance ...\n');
end
if wts_passed == 0
    %
    %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    meanvar  = zeros(ncomps,1);      % size of the projections
    if ncomps == urchans % if weights are square . . .
        winv = inv(weights*sphere);
    else
        fprintf('Using pseudo-inverse of weight matrix to rank order component projections.\n');
        winv = pinv(weights*sphere);
    end
    for s=1:ncomps
        if verbose,
            fprintf('%d ',s);         % construct single-component data matrix
        end
        % project to scalp, then add row means
        compproj = winv(:,s)*activations(s,:);
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
    windex = windex(ncomps:-1:1); % order large to small
    meanvar = meanvar(windex);
    %
    %%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
    %
    if nargout>2, % if activations are to be returned
        if verbose,
            fprintf('Permuting the activation wave forms ...\n');
        end
        activations = activations(windex,:);
    else
        clear activations
    end
    weights = weights(windex,:);% reorder the weight matrix
    bias  = bias(windex);		% reorder them
    signs = diag(signs);        % vectorize the signs matrix
    signs = signs(windex);      % reorder them
else
    fprintf('Components not ordered by variance.\n');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

return

if nargout > 6
    u=weights*data + bias*ones(1,frames);
    y = zeros(size(u));
    for c=1:chans
        for f=1:frames
            y(c,f) = 1/(1+exp(-u(c,f)));
        end
    end
end
