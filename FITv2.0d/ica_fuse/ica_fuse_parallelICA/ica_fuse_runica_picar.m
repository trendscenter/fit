%
%runica_parallelicaMul() - Perform Parallal Independent Component Analysis
%   ICA decomposition,  based on  informax ICA decomposition code.
%   the infomax ICA algorithm of Bell & Sejnowski (1995) with the natural gradient feature
%   of Amari, Cichocki & Yang, the extended-ICA algorithm of Lee, Girolami & Sejnowski, PCA dimension reduction,
%   and/or specgram() preprocessing (M. Zibulevsky).
%
%Parallel ICA is one type of constrained ICA, with the correlation between
%  two data emphasized. In this particular code the correlation is defined as
%  correlation between multiple columns of A matrix from data 1 and columns of A matrix from 
%  data 2. The columes are selected by surviving a thredhold during the optimazation process
% The square of the correlation is the constrain term.
% 
%
% runica_parallelicaMul() - Perform Parallel Independent Component Analysis (ICA) decomposition
%            of two psychophysiological datasets using the infomax ICA algorithm of
%            Bell & Sejnowski (1995) with the natural gradient feature
%            of Amari, Cichocki & Yang, the extended-ICA algorithm
%            of Lee, Girolami & Sejnowski, PCA dimension reduction,
%            and/or specgram() preprocessing (M. Zibulevsky).
% Usage:
%      simply >> [weights,sphere] = runica_parallelicaMul(data,'dewhitem',{dewhiteMfMRI},'whitem',{whiteMfMRI},'constrained_components',1);
%       or
%        else >> [weights,sphere,activations,bias,signs,lrates] ...
%                                 = runica_parallelicaMul(data,'Key1',Value1,...);
% Input_Variable:
%
% data        = cell type input data:{data1, data2};data1:(chans,frames*epochs);data2:(chans,frames*epochs).
%               Note: If data consists of multiple discontinuous epochs,
%               each epoch should be separately baseline-zero'd using:
%                  >> data = rmbase(data,frames,basevector);
% Key Keywods value
% 'dewhitem'    = cell type Dewhite matrix {dewhiteM1, dewhiteM2}.. 
% 'whitem'      = cell type white matrix {whiteM1, whiteM2}.. 
% 'constrained_components'  = the maximum pair of components being constrained 
% 'constrained_connection'  = the threshold, the connection above which would be strengthed
% 'endurance'   = number ,the maximumlly allowed descend trend of entropy 
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
%               every N training blocks. If N < 0, fix number of
%               sub-Gaussian
%               components to -N [faster than N>0]      (default|0 -> off)
% 'specgram'  = [srate loHz hiHz Hzinc frames] decompose a complex time/frequency
%               transform of the data         (defaults [srate 0 srate/2 1 0])
% 'posact'    = make all component activations net-positive(default 'on'}
% 'verbose'   = give ascii messages ('on'/'off')        (default -> 'on')
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
%  10-07-97 put acos() outside verbose loop; verbose 'off' wasn't stopping -sm
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

function [weights,sphere,activations,tf_reorg,c,indep,dist_r, bias,signs,lrates,y] = ica_fuse_runica_picar(data,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)

if nargin < 5
    help ica_fuse_runica_picar;
    return
end

% separte two datasets from input argument
data1=data{1};
data2=data{2};
clear data;

[chans(1) frames(1)] = size(data1); % determine the data size
urchans(1) = chans(1);  % remember original data channels
datalength(1) = frames(1);

[chans(2) frames(2)] = size(data2); % determine the data size
urchans(2) = chans(2);  % remember original data channels
datalength(2) = frames(2);


%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 3;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.9;      % anneal by multiplying lrate by this original 0,9 to 0.95 changed by JL
DEFAULT_EXTANNEAL    = 0.98;      % or this if extended-ICA
DEFAULT_MAXSTEPS     = 512;       % ]top training after this many steps 512
DEFAULT_MOMENTUM     = 0.0;       % default momentum weight

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.9;       % if weights blowup, restart with lrate
% lower by this factor
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
DEFAULT_LRATE        = 0.015./log(chans);

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

DEFAULT_SPHEREFLAG   = 'on';      % use the sphere matrix as the default
%   starting weight matrix
DEFAULT_PCAFLAG      = 'off';     % don't use PCA reduction
DEFAULT_POSACTFLAG   = 'on';      % use posact()
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule

%--constrained ICA parameters
CONSTRAINED_COMPONENTS=3; % NUMBER OF COMPONENTS FROM EACH DATASET BEING CONSTRAINED
CONSTRAINED_CONNECTION=0.5; % CORRELATION THRESHOLD TO BE CONSTRAINED; HIGH THRESHOLD WILL BE STRENGTHENED.
ENDURANCE = -1e-3; % the maximumlly allowed descending trend of entropy;

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
kurtsize   = [MAX_KURTSIZE,MAX_KURTSIZE];
signsbias  = 0.02;                   % bias towards super-Gaussian components
extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
nsub       = DEFAULT_NSUB;
wts_blowup = [0,0];                      % flag =1 when weights too large
wts_passed = 1;                      % flag weights passed as argument
%                 
Connect_threshold =CONSTRAINED_CONNECTION; % set a threshold to select columns constrained.
MaxComCon  =       CONSTRAINED_COMPONENTS;
trendPara  = ENDURANCE; %depends on the requirement on connection; the more negative,the stronger the contrains ,that may cause overfitting

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
        MaxComCon  = Value ;
    elseif strcmp(Keyword,'constrained_connection')        
        Connect_threshold =Value; % set a threshold to select columns constrained.                
    elseif strcmp(Keyword,'endurance')
        trendPara =Value; %
    elseif strcmp(Keyword,'ref_snp');
        ref = Value;      
    else
        fprintf('runica(): unknown flag')
        return
    end
end
%
% 
% if ~exist('dewhiteM','var') | ~exist('whiteM','var')
%     help ica_fuse_runica_parallelica2
%     return
% end

% if ~isempty(ica_options)
%     Connect_threshold = ica_options{4}; % set a threshold to select columns constrained.
%     MaxComCon = ica_options{2};
%     trendPara  = ica_options{6}; %depends on the requirement on connection; the more negative,the stronger the contrains ,that may cause overfitting
% end

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
%if verbose,
%    fprintf('Removing mean of each channel ...\n');
%end
disp('Not removing mean of each channel!!!');
%data = data - mean(data')'*ones(1,frames);      % subtract row means

if verbose,
    fprintf('Final training data1 range: %g to %g\n', ...
        min(min(data1)),max(max(data1)));
    fprintf('Final training data2 range: %g to %g\n', ...
        min(min(data2)),max(max(data2)));
    
end
%
%%%%%%%%%%%%%%%%%%% Perform PCA reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Reducing the data to %d principal dimensions...\n',ncomps);
    [eigenvectors1,eigenvalues1,data1] = pcsquash(data1,ncomps(1)); % changed for two dsatasets
    [eigenvectors2,eigenvalues2,data2] = pcsquash(data2,ncomps(2)); % changed for two datasets
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
    sphere{1} = 2.0*inv(sqrtm(cov(data1'))); % find the "sphering" matrix = spher()
    sphere{2} = 2.0*inv(sqrtm(cov(data2'))); % find the "sphering" matrix = spher()
    if ~weights,
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
        end
        weights=[];
        weights{1} = eye(ncomps(1),chans(1)); % begin with the identity matrix
        weights{2} = eye(ncomps(2),chans(2)); % begin with the identity matrix
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
elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~weights
        if verbose,
            fprintf('Using the sphering matrix as the starting weight matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere{1} = 2.0*inv(sqrtm(cov(data1'))); % find the "sphering" matrix = spher()
        weights{1} = eye(ncomps(1),chans(1))*sphere{1}; % begin with the identity matrix
        sphere{1} = eye(chans(1));                 % return the identity matrix
        sphere{2} = 2.0*inv(sqrtm(cov(data2'))); % find the "sphering" matrix = spher()
        weights{2} = eye(ncomps(2),chans(2))*sphere{2}; % begin with the identity matrix
        sphere{2} = eye(chans(2));                 % return the identity matrix
        
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere{1} = eye(chans(1));                 % return the identity matrix
        sphere{2} = eye(chans(2));                 % return the identity matrix
    end
elseif strcmp(sphering,'none')
    sphere{1} = eye(chans(1));                     % return the identity matrix
    sphere{2} = eye(chans(2));                     % return the identity matrix
    if ~weights
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        weights{1} = eye(ncomps(1),chans(1)); % begin with the identity matrix
        weights{2} = eye(ncomps(2),chans(2)); % begin with the identity matrix
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
    end
    sphere{1} = eye(chans(1),chans(1));
    sphere{2} = eye(chans(2),chans(2));
    if verbose,
        fprintf('Returned variable "sphere" will be the identity matrix.\n');
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%
lastt=fix((datalength./block-1).*block+1);
BI{1}=block(1)*eye(ncomps(1),ncomps(1));
BI{2}=block(2)*eye(ncomps(2),ncomps(2));
delta{1}=zeros(1,chans(1)*ncomps(1));
delta{2}=zeros(1,chans(2)*ncomps(2));
changes = [];
degconst = 180./pi;
startweights = weights;
prevweights = startweights;
oldweights = startweights;
prevwtchange{1} = zeros(chans(1),ncomps(1));
prevwtchange{2} = zeros(chans(2),ncomps(2));
oldwtchange{1} = zeros(chans(1),ncomps(1));
oldwtchange{2} = zeros(chans(2),ncomps(2));

lrates = zeros(2,maxsteps);
onesrow{1} = ones(1,block(1));
onesrow{2} = ones(1,block(2));
bias{1} = zeros(ncomps(1),1);
bias{2} = zeros(ncomps(2),1);
signs{1} = ones(1,ncomps(1));    % initialize signs to nsub -1, rest +1
signs{2} = ones(1,ncomps(2));    % initialize signs to nsub -1, rest +1
for k=1:nsub
    signs{1}(k) = -1;
    signs{2}(k) = -1;
end
if extended & extblocks < 0 & verbose,
    fprintf('Fixed extended-ICA sign assignments:  ');
    for k=1:ncomps
        fprintf('%d ',signs(k));
    end; fprintf('\n');
end
signs{1} = diag(signs{1}); % make a diagonal matrix
signs{2} = diag(signs{2}); % make a diagonal matrix

oldsigns{1} = zeros(size(signs{1}));;
oldsigns{2} = zeros(size(signs{2}));;
change=[0,0];
signcount =[ 0,0];              % counter for same-signs
signcounts = [];
urextblocks = extblocks;    % original value, for resets
old_kk{1} = zeros(1,ncomps(1));   % for kurtosis momemtum
old_kk{2} = zeros(1,ncomps(2));   % for kurtosis momemtum

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
step=[0,0];
blockno =[1,1];  % running block counter for kurtosis interrupts
stopsign=[0,0];      angledelta=[0,0];
alphak=[1,1];%alphak_R=[1,1];
rrate=[1,1];

data_t = data2';
taget_index = [];
null_index = [];
clear dist_r;
for j1 = 1:size(ref,1)
    target_index{j1} = find(ref(j1,:));
    null_index{j1} = find(~ref(j1,:));
    ref(j1,:) = ref(j1,:)/norm(ref(j1,:));
    dist_r{j1} = [];
end
% temp1 = zeros(1,length(ref));
% temp1(1:1000) = 1;
% temp1 = temp1/norm(temp1);
% ref(target_index) = temp1(1);

indep = [];
weight_mo1 = 1;
weight_mo2 = 9*ones(1,size(ref,1));
tf_ref = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while step(1) < maxsteps | step(2) < maxsteps
    
    %  data1 ICA training 
    if (~stopsign(1)) 
        
        % the beginning entropy of each step
        u=weights{1}*data1(:, :) + bias{1}*ones(1,frames(1));
        y=1./(1+exp(-u));
        temp=log(abs(weights{1}*y.*(1-y))+eps); entropy(1)=mean(temp(:));
        %------------
        
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
        
        for j1 = 1:size(weights{1},1)
            weights{1}(j1,:) = weights{1}(j1,:)/norm(weights{1}(j1,:));
        end
       
        %---------------------------
        % if weight is not  blowup, update
        if ~wts_blowup(1)
            
            oldwtchange{1} = weights{1}-oldweights{1};
            step(1)=step(1)+1;
            lrates(1,step(1)) = lrate(1);
            angledelta(1)=[0];        
            delta{1}=reshape(oldwtchange{1},1,chans(1)*ncomps(1));
            change(1)=delta{1}*delta{1}'       ;
            
        end
        %DATA1 blow up restart-------------------------------        
        if wts_blowup(1) | isnan(change(1))|isinf(change(1)),  % if weights blow up,
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
                oldsigns{1} = zeros(size(signs{1}));;
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
            u=weights{1}*data1(:, :) + bias{1}*ones(1,frames(1));
            y=1./(1+exp(-u));
            temp=log(abs(weights{1}*y.*(1-y))+eps);
            lossf1(step(1))=mean(temp(:));
            %changes of entropy term added by jingyu
            if step(1)>1
                entropychange(1)=lossf1(step(1))-entropy(1);
            else
                entropychange(1)=1;
            end
            %--------
            
            if step(1)>5 
                index1 = ica_fuse_falsemaxdetect(lossf1,trendPara); 
                if index1
                    rrate(1) = rrate(1)*0.9;         % anneal learning rate empirical
                end
            end % end of test------------------------
            
            
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step(1)> 2 & ~stopsign(1)
                angledelta(1)=acos((delta{1}*olddelta{1}')/sqrt(change(1)*oldchange(1)));
            end
            if verbose,
                if step(1) > 2,
                    if ~extended,
                        fprintf(...
                            'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\t\t', ...
                            step(1),lrate(1),change(1),degconst*angledelta(1));
                    else
                        fprintf(...
                            'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\t\t',...
                            step(1),lrate(1),change(1),degconst*angledelta(1),(ncomps(1)-sum(diag(signs{1})))/2);
                    end
                elseif ~extended
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f\t\t',step(1),lrate(1),change(1));
                else
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f, %d subgauss\t\t',...
                        step(1),lrate(1),change(1),(ncomps(1)-sum(diag(signs{1})))/2);
                end % step > 2
            end; % if verbose
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if degconst*angledelta(1) > annealdeg  % entropychange(1)<0  %  % 
                lrate(1) = lrate(1)*annealstep;          % anneal learning rate
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
            
        end; % end if weights in bounds
        
    end
    
    
    %----------------
    %  data2 ICA training 
    if ~stopsign(2)
        
        % the beginning entropy of each step
        u=weights{2}*data2(:, :) + bias{2}*ones(1,frames(2));
        y=1./(1+exp(-u));
        temp=log(abs(weights{2}*y.*(1-y))+eps); entropy(2)=mean(temp(:));
        %------------
        
        % begin to update W matrix        
        permuteVec = randperm(datalength(2)); % shuffle data order at each step
        for t=1:block(2):lastt(2), %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
            data2_temp = data2(:, permuteVec(t:t+block(2)-1));
            if biasflag
                
                u=weights{2}*data2_temp + bias{2}*onesrow{2};
                
            else
                u=weights{2}*data2_temp;
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
        
        temp1 = weights{2};
        for j1 = 1:size(weights{2},1)
            temp1(j1,:) = temp1(j1,:)/norm(temp1(j1,:));
        end
        temp1_p = temp1 - oldweights{2};
        temp2 = [];
        for j1 = 1:size(weights{2},1)
            temp2(j1,1) = norm(temp1_p(j1,:));
        end
        
       
        icasig_temp = oldweights{2}*data2;
        for j1 = 1:size(weights{2},1)
            icasig_temp(j1,:) = icasig_temp(j1,:)/norm(icasig_temp(j1,:));
        end

        
        if step(2) >= 1
            
            weights_orig = weights;
            comp_index_all = [];
            for k1 = 1:size(ref,1)
                y1 = abs(icasig_temp(:,target_index{k1})) - repmat(ref(k1,target_index{k1}),size(weights{2},1),1);
                d1 = sum(y1.^2,2);
                [mindist,comp_index] = sort(d1);
                comp_index_all(k1,:) = comp_index(1:2);
                
                sx = (weights_orig{2}' \ dewhiteM{2}')';

                weights_temp = weights_orig;
                j1 = comp_index(1);
                fk = sum(y1(j1,:).^2);
                deltafk = (abs(icasig_temp(j1,target_index{k1}))-ref(k1,target_index{k1}))*(repmat(sign(icasig_temp(j1,target_index{k1})'),1,size(weights{2},2)).*data_t(target_index{k1},:))*temp1'*temp1;
                P = deltafk/norm(deltafk);
                weights{2}(j1,:) = weights{2}(j1,:) - weight_mo2(k1)*temp2(j1)*P;
                weights_temp{2}(j1,:) = weights_temp{2}(j1,:) - weight_mo2(k1)*temp2(j1)*P;
                
                ref_index(k1) = j1;
                
                icasig_temp1 = weights_temp{2}*data2;
                icasig_temp1 = icasig_temp1(ref_index(k1),:)/norm(icasig_temp1(ref_index(k1),:));
                phi2 = sum(sum(abs(weights_temp{2}*data2)))/size(weights{2},1);
                phi3 = sum((abs(icasig_temp1(target_index{k1})) - ref(k1,target_index{k1})).^2);
                dist_r{k1} = [dist_r{k1};phi3];
                
                indep_temp = log(abs(det(weights_temp{2})));
                if step(2) > 10
                    if indep_temp < indep(end)
                        weight_mo2(k1) = 0.95*weight_mo2(k1);
                    elseif dist_r{k1}(end) > dist_r{k1}(end-1)
                        weight_mo2(k1) = min(9,1.01*weight_mo2(k1));
                    else
                        temp1 = (dist_r{k1}(end-1)^2)/(dist_r{k1}(end)*dist_r{k1}(end-2));
                        weight_mo2(k1) = weight_mo2(k1)*temp1;
                    end
                end
                
            end
            
            % normalization
            
            for j1 = 1:size(weights{2},1)
                weights{2}(j1,:) = weights{2}(j1,:)/norm(weights{2}(j1,:));
            end
            
            phi1 = log(abs(det(weights{2})));
            indep = [indep;phi1];
        end
        
        if ~tf_ref
            if step(2) == 10
                tf_ref = 1;
                c = picamr_cluster(comp_index_all);
                tf_reorg = 0;
                for j1 = 1:length(c)
                    if length(c{j1}) > 1
                        tf_reorg = 1;
                        break;
                    end
                end
                
                if tf_reorg
                    display('reorganize reference snps');
                    wts_blowup(2) = 1;
                    % wts_blowup(1) = 1;
                    ref_old = ref;
                    ref = zeros(length(c),size(ref_old,2));
                    for j1 = 1:length(c)
                        for j2 = 1:length(c{j1})
                            ref(j1,find(ref_old(c{j1}(j2),:))) = 1;
                        end
                    end
                    taget_index = [];
                    null_index = [];
                    clear dist_r;
                    for j1 = 1:size(ref,1)
                        target_index{j1} = find(ref(j1,:));
                        null_index{j1} = find(~ref(j1,:));
                        ref(j1,:) = ref(j1,:)/norm(ref(j1,:));
                        dist_r{j1} = [];
                    end
                end
            end
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
            u=weights{2}*data2(:, :) + bias{2}*ones(1,frames(2));
            y=1./(1+exp(-u));
            temp=log(abs(weights{2}*y.*(1-y))+eps);
            lossf2(step(2))=mean(temp(:));
            %changes of entropy term added by jingyu
            if step(2)>1;
                entropychange(2)=lossf2(step(2))-entropy(2);
            else
                entropychange(2)=1;
            end
            %--------
            if step(2)>5 
                index1 = ica_fuse_falsemaxdetect(lossf2,trendPara); 
                if index1
                    rrate(2) = rrate(2)*0.9;         % anneal learning rate empirical
                end
            end % end of test------------------------
            
            
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step(2)> 2 & ~stopsign(2)
                angledelta(2)=acos((delta{2}*olddelta{2}')/sqrt(change(2)*oldchange(2)));
            end
            
            
            if verbose,
                if step(2) > 2,
                    if ~extended,
                        fprintf(...
                            'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
                            step(2),lrate(2),change(2),degconst*angledelta(2));
                    else
                        fprintf(...
                            'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\n',...
                            step(2),lrate(2),change(2),degconst*angledelta(2),(ncomps(2)-sum(diag(signs{2})))/2);
                    end
                elseif ~extended
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f\n',step(2),lrate(2),change(2));
                else
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f, %d subgauss\n',...
                        step(2),lrate(2),change(2),(ncomps(2)-sum(diag(signs{2})))/2);
                end % step > 2
            end; % if verbose
            
            
            %           
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if degconst*angledelta(2) > annealdeg % entropychange(2)<0 
                lrate(2) = lrate(2)*annealstep;          % anneal learning rate
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
            
        end; % end if weights in bounds        
    end
    
    % ------------------------- 
    % modifying weights based on correlation between data1 A Matrix and data2 A Matrix 
    % %%%%%%%%%%%%%%%%%nudging
    
    if (step(1)>10 & step(2)>10) & ( ~stopsign(1) | ~stopsign(2)) 
        
        % A matrix of data1 
        mx = (weights{1}' \ dewhiteM{1}')';
        %  A matrix of data2, 
        sx = (weights{2}' \ dewhiteM{2}')';
        % calculate the correlation of all componentsts   % match tow modality components
        for k1 = 1:size(ref,1)
            [Corr_matrix]=abs(ica_fuse_corr(mx,sx));
            maxcorr=[];maxcol=[];maxrow=[];
            [maxcol1,maxrow1] = find(Corr_matrix == max(Corr_matrix(:,ref_index(k1))));
            maxcorr1 = Corr_matrix(maxcol1,maxrow1);
            if MaxComCon == 1
                if maxcorr1 > Connect_threshold;
                    ix = 0;
                    maxcol = [maxcol1];
                    maxrow = [maxrow1];
                    maxcorr = [maxcorr1];
                else
                    ix = [];
                    maxcorr = maxcorr1;
                end
            else
                Corr_matrix(:,ref_index(k1)) = 0;
                for j=1:min([ncomps(1),ncomps(2)])
                    [mm,ss]=find(Corr_matrix==max(Corr_matrix(:)));
                    if length(mm) == 1
                        maxcol(j)=mm; maxrow(j)=ss; maxcorr(j)=Corr_matrix(mm,ss);
                        Corr_matrix(mm,:)=0;Corr_matrix(:,ss)=0;
                    end
                end
                %[temp,index]=sort(abs(maxcorr),'descend');
                ix=find(maxcorr>Connect_threshold);
                if length(ix)>MaxComCon-1; ix=ix(1:MaxComCon-1); end
                if maxcorr1 > Connect_threshold;
                    ix = [0,ix];
                    maxcol = [maxcol1,maxcol];
                    maxrow = [maxrow1,maxrow];
                    maxcorr = [maxcorr1,maxcorr];
                end
            end
        end
                   
        if ~ isempty(ix)
            
            for Cons_com=1:length(ix); % constraned componenets
                
                
                % Updata the weights
                a=mx(:,maxcol(Cons_com))';u=sx(:,maxrow(Cons_com))'; % 
                b=cov(a,u);  b=b(2); tmcorr=b/std(a)/std(u);
                comterm=2*b/var(a)/var(u); coml=length(a);
                
                if ~stopsign(1) %& ~Overindex1
                    
                    deltaf=comterm*(u-mean(u)+b*(mean(a)-a)/var(a));% 1st order derivative
                    P=deltaf./norm(deltaf); %;(H1*deltaf')';                          
                    alphak(1)=findsteplength3(-tmcorr^2,-deltaf,a,u,alphak(1),P,0.0001,0.999);
                    aweights_temp=1e3*lrate(1)*rrate(1)*alphak(1)*P; % 
                    mx(:,maxcol(Cons_com))=mx(:,maxcol(Cons_com)) + aweights_temp';
                end
                
                if ~stopsign(2) % & ~Overindex2 
                    deltaf=(comterm*( a-mean(a) + b/var(u)*(mean(u) -u)));% 1st order derivative
                    P=deltaf./norm(deltaf) ;% (H2*deltaf')';                          
                    alphak(2)=findsteplength3(-tmcorr^2,-deltaf,u,a,alphak(2),P,0.0001,0.999);
                    aweights_temp=1e3*lrate(2)*rrate(2)*alphak(2)*P; % 
                    sx(:,maxrow(Cons_com))=sx(:,maxrow(Cons_com)) + aweights_temp';
                    
                end
            end
            
            if ~stopsign(1)
                temp=weights{1};
                weights{1} = whiteM{1}*mx;   weights{1} =inv(weights{1});% weights{1} \ eye(size(weights{1}));
                if max(max(abs(weights{1}))) > MAX_WEIGHT
                    rrate(1)=rrate(1)*0.95;
                    weights{1}=temp;
                end
            end
            
            if  ~stopsign(2)
                temp=weights{2};
                weights{2} = whiteM{2}*sx;
                weights{2} =inv(weights{2});% weights{1} \ eye(size(weights{1}));
                if max(max(abs(weights{2}))) > MAX_WEIGHT
                    rrate(2)=rrate(2)*0.95;
                    weights{2}=temp;
                end
            end
            %             test ---------------------
            
            sx = (weights{2}' \ dewhiteM{2}')';
            mx = (weights{1}' \ dewhiteM{1}')';
%             sx = dewhiteM{2}*pinv(weights{2}*sphere{2});
%             mx = dewhiteM{1}*pinv(weights{1}*sphere{1});            
            a=mx(:,maxcol(1)); b=sx(:,maxrow(1));
            temp=corrcoef(a,b);  temp=temp(2);% correlation
            
            % normalization
            for j1 = 1:size(weights{1},1)
                weights{1}(j1,:) = weights{1}(j1,:)/norm(weights{1}(j1,:));
            end
            for j1 = 1:size(weights{2},1)
                weights{2}(j1,:) = weights{2}(j1,:)/norm(weights{2}(j1,:));
            end
        
            if abs(temp)<maxcorr(1);
                disp('Wrong direction !!!! ');    
            end
            lossf3(max(step))=abs(temp);                   
            %             -----------------end test
            oldweights{2} = weights{2};
            oldweights{1} = weights{1};
        end
        
    end%--------------------------
    
    if stopsign(1)==1 & stopsign(2)==1
        laststep=step;
        step=[maxsteps,maxsteps];                % stop when weights stabilize
    end
    
    
end; % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('laststep', 'var')
    lrates = lrates(:,1:max(laststep));           % truncate lrate history vector
end

%figure;plot(lossf1);
%figure;plot(lossf2);
%if exist('lossf3','var')
%    figure;plot(lossf3);
%end


%%%%%%%%%%%%%% Orient components towards positive activation %%%%%%%%%%%
%
if strcmp(posactflag,'on')
    [activations{1},winvout{1},weights{1}] = ica_fuse_posact(data1,weights{1});
    [activations{2},winvout{2},weights{2}] = ica_fuse_posact(data2,weights{2});
    % changes signs of activations and weights to make activations
    % net rms-positive
else
    activations{1} = weights{1}*data1;
    activations{2} = weights{2}*data2;
end
%
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
end
%
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
        winv = inv(weights{1}*sphere{1});
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
        winv = inv(weights{2}*sphere{2});
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

function gamma = findsteplength1(fk,deltafk,x1,x2,gamma,P,c1,c2,lrate,null_index)
% 0<c1<0.5<c2<1;
% (f(x+ap)-f(x))/c1/dela_f(x)/P <a
con=1;coml=length(x1);
while con & con<100
    xnew=x1-lrate*gamma*deltafk;   % new weights 
    for j1 = 1:size(xnew,1)
        xnew(j1,:) = xnew(j1,:)/norm(xnew(j1,:));
    end
    y1 = xnew*x2;    % new components
    for j1 = 1:size(xnew,1)
        y1(j1,:) = y1(j1,:)/norm(y1(j1,:));
    end
    fk1=sum(sum(abs(y1(:,null_index)))); 
    deltafk1=sign(y1(:,null_index))*x2(:,null_index)';% 1st order derivative
    P_t = P';
    firstterm1=(fk1-fk)/c1;firstterm2=trace(deltafk*P_t);
    secondterm1=trace(deltafk1*P_t);secondterm2=trace(deltafk*P_t);
    
    if firstterm1 > lrate*gamma*firstterm2
        if firstterm1<0;
            gamma = 0.9*abs(firstterm1/firstterm2); 
        else
            gamma = 0.5*gamma;
         
        end
        if gamma<1e-6 
            gamma=0;con=0;
        else
            con=con+1;
        end
        
    elseif secondterm1<0 & secondterm1<c2*secondterm2  
        
        % alphak=abs(secondterm1/secondterm2/c2)*alphak;
        con=con+1;
        gamma=1.1*gamma;     
    elseif secondterm1 >0 & secondterm2 <0  
        gamma=0.9*gamma;    con=con+1;
        
    else
        con=0;
    end
end
if con>=50 gamma=0; disp('learning rate searching for L1 norm failed!'); end

function gamma = findsteplength2(fk,deltafk,x1,x2,x3,gamma,P,c1,c2,lrate,target_index)
% 0<c1<0.5<c2<1;
% (f(x+ap)-f(x))/c1/dela_f(x)/P <a
con=1;coml=length(x1);
old_firstterm1 = 1e6;
a1 = 1.05;
while con & con<100
    xnew = x1 - lrate*gamma*deltafk;   % new weights 
    xnew = xnew/norm(xnew);
    y1 = xnew*x3;    % new components
    y1 = y1/norm(y1);
    fk1 = sum((abs(y1(target_index))-x2(target_index)).^2); 
    deltafk1 = (abs(y1(target_index))-x2(target_index))*...
        (repmat(sign(y1(target_index)'),1,size(x1,2)).*x3(:,target_index)'); % 1st order derivative
    P_t = P';
    firstterm1=(fk1-fk)/c1;firstterm2=trace(deltafk*P_t);
    secondterm1=trace(deltafk1*P_t);secondterm2=trace(deltafk*P_t);
    
    if firstterm1 > lrate*gamma*firstterm2
        if firstterm1<0;
            gamma = 0.9*abs(firstterm1/firstterm2); 
        else
            if firstterm1 < old_firstterm1
                gamma = a1*gamma;
            else
                gamma = (1/(a1+0.1))*gamma;
            end
            old_firstterm1 = firstterm1;
        end
        if gamma<1e-6 
            gamma=0;con=0;
        else
            con=con+1;
        end
        
    elseif secondterm1<0 & secondterm1<c2*secondterm2  
        
        % alphak=abs(secondterm1/secondterm2/c2)*alphak;
        con=con+1;
        gamma=1.1*gamma;     
    elseif secondterm1 >0 & secondterm2 <0  
        gamma=0.9*gamma;    con=con+1;
        
    else
        con=0;
    end
end
if con>=50 gamma=0; disp('learning rate searching for dist failed!'); end

function alphak=findsteplength3(fk,deltafk,x1,x2,alphak,P,c1,c2)
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
if con>=50 alphak=0; disp('learning rate searching for fMRI data failed!'); end


% clustering of dynamically selected constrained components for pica-r

function c = picamr_cluster(comp_index)

c = [];
k = 1;
ind = 1:size(comp_index,1);
while ~isempty(comp_index)
    r = comp_index(1,:);
    ind_comb1 = find(comp_index(:,1) == r(1));
    ind_comb2 = find((comp_index(:,2) == r(1)) | (comp_index(:,1) == r(2)));
    ind_comb = union(ind_comb1,ind_comb2);
    c{k} = ind(ind_comb);
    comp_index(ind_comb,:) = [];
    k = k+1;
    ind(ind_comb) = [];
end

    
