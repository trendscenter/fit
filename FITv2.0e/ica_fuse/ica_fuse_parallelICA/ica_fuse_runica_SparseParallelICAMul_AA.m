%
%ica_fuse_runica_SparseParallelICAMul_AA() - Perform Sparse Parallal Independent Component Analysis decomposition based on informax ICA and Hoyer projection.
%
%Sparse Parallel ICA is one type of constrained ICA, with the correlation between two data emphasized and sparsity of sources being enhanced. 
%In this particular code the correlation is defined as correlation between multiple columns of A matrix from data 1 and columns of A matrix from data 2. 
%The columes are selected if passing a correlation thredhold during the optimazation process
%The square of the correlation is the constrain term.
%The sparsity of a source is optimized when its hoyer index is smaller than the preset Hoyer value

% Usage:
%      simply >> [weights,sphere,activations] = ica_fuse_runica_SparseParallelICAMul_AA(data,'dewhitem', dewhiteM, 'whitem', whiteM,'constrained_components',1,'hoyerindex',[0.1,0.45]);
%       or
%        else >> [weights,sphere,activations] = ica_fuse_runica_SparseParallelICAMul_AA(data,'Key1',Value1,...);
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
function [weights,sphere_t,activations,laststep,sparse_sign,max_corr_corropt,bias,signs,lrates,y] = ica_fuse_runica_SparseParallelICAMul_AA(data,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)

if nargin < 5
    help ica_fuse_runica_SparseParallelICAMul_AA
    return
end

% separte two datasets from input argument
data1=data{1};
data2=data{2};
data1_pca = data1;
data2_pca = data2;
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
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      % anneal by multiplying lrate by this original 0,9 to 0.95 changed by JL
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

DEFAULT_SPHEREFLAG   = 'none';%% 'on';       % use the sphere matrix as the default
%   starting weight matrix
DEFAULT_PCAFLAG      = 'off';     % don't use PCA reduction
DEFAULT_POSACTFLAG   = 'on';      % use posact()
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule

%--constrained ICA parameters
CONSTRAINED_COMPONENTS = 1; % NUMBER OF COMPONENTS FROM EACH DATASET BEING CONSTRAINED
CONSTRAINED_CONNECTION = 0.18; % CORRELATION THRESHOLD TO BE CONSTRAINED; HIGH THRESHOLD WILL BE STRENGTHENED.
ENDURANCE = -1e-5; % the maximumlly allowed descending trend of entropy;

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
wts_passed = 0;                      % flag weights passed as argument
%                 
Connect_threshold = CONSTRAINED_CONNECTION; % set a threshold to select columns constrained.
MaxComCon  =       CONSTRAINED_COMPONENTS;
trendPara  = ENDURANCE; %depends on the requirement on connection; the more negative,the stronger the contrains,that may cause overfitting
%
z_thr = 3.5;
SBR_min_sparstywork = 3;
SBR_min = [5,5];    %%minimum SBR value
delta_lambda_hoyer(1) = 0.005;  %%step size for hoyer value for modality 1
delta_lambda_hoyer(2) = 0.0005; %%step size for hoyer value for modality 2
sparsty_opt_start_step = 8; %%the step to start sparsity optimization
corr_opt_start_step = 5; %%the step to start correlation optimization
sparse_sign1_clct = [];
sparse_sign2_clct = [];

mx_degree = [];
sx_degree = [];

corr_angle_thr = 30;%%threshold for the angle between a component loading before correlation optimization and a component loading after correlation optimization
change_flag_mx = [];
change_flag_sx = [];

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
    elseif strcmp(Keyword,'hoyerindex')
      if isstr(Value)
         fprintf('runica(): Hoyer index value must be between 0 and 0.9')
         return
      end
      if Value(1)<0 || Value(1) >0.9 || Value(2)<0 || Value(2) >0.9
         fprintf('runica(): Hoyer index value must be between 0 and 0.9');
         return
      end
      lambda_hoyer1 = Value(1);
      lambda_hoyer2 = Value(2);
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
%     sphere{1} = 2.0*inv(sqrtm(cov(data1'))); % find the "sphering" matrix = spher()
%     sphere{2} = 2.0*inv(sqrtm(cov(data2'))); % find the "sphering" matrix = spher()
    [v{1}]=eig(data1*data1'); sphere{1}=1./mean(sqrt(v{1}))*eye(size(data1,1)); % need pre-whiten done
    [v{2}]=eig(data2*data2'); sphere{2}=1./mean(sqrt(v{2}))*eye(size(data2,1)); % need pre-whiten done
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
        weights = [];
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
Crate=[1,1];
corr_hi_track = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corropt_step = 1;
change_flag1 = zeros(maxsteps,1);
change_flag2 = zeros(maxsteps,1);
sphere_t{1} = sphere{1};
sphere_t{2} = sphere{2};
max_corr_sparseopt2 = [];
max_corr_corropt = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while step(1) < maxsteps | step(2) < maxsteps
    
    %  data1 ICA training 
    if ~stopsign(1)
        
        % the beginning entropy of each step
        entropy(1) = Compute_entropy(weights{1},data1);

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
            change_flag1 = zeros(maxsteps,1);
            corropt_step = 1;
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
            %test the trend of entropy term, avoid the overfitting of correlation
            lossf1(step(1)) = Compute_entropy(weights{1},data1);
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
                    Crate(1) = Crate(1)*0.9;    %%anneal learning rate empirical
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
                            'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
                            step(1),lrate(1),change(1),degconst*angledelta(1));
                    else
                        fprintf(...
                            'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\n',...
                            step(1),lrate(1),change(1),degconst*angledelta(1),(ncomps(1)-sum(diag(signs{1})))/2);
                    end
                elseif ~extended
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f\n',step(1),lrate(1),change(1));
                else
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f, %d subgauss\n',...
                        step(1),lrate(1),change(1),(ncomps(1)-sum(diag(signs{1})))/2);
                end % step > 2
            end; % if verbose
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if entropychange(1)<0 %|degconst*angledelta(1) > annealdeg, %% 
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
        entropy(2) = Compute_entropy(weights{2},data2);
        %------------
        
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
            change_flag2 = zeros(maxsteps,1);
            corropt_step = 1;
            if extended
                signs{2} = ones(1,ncomps(2));    % initialize signs to nsub -1, rest +1
                for k=1:nsub
                    signs{2}(k) = -1;
                end
                signs{2} = diag(signs{2}); % make a diagonal matrix
                oldsigns{2} = zeros(size(signs{2}));;
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
            
            %test the trend of entropy term, avoid the overfitting of correlation
            lossf2(step(2)) = Compute_entropy(weights{2},data2);
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
                    Crate(2) = Crate(2)*0.9;  %%anneal learning rate empirical
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
            if entropychange(2)<0 %|degconst*angledelta(2) > annealdeg, %%
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

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%----------------------Sparsity optimization----------------------
    %%%%Perform Hoyer projection on two modilaties to optimize sparsity and then update the data
    data1_org = data1;
    data2_org = data2;
    icasig1 = weights{1}*data1;
    icasig2 = weights{2}*data2;
    
    mx_tmp = inv(weights{1});
    mx = (dewhiteM{1}*mx_tmp);
    %%%  A matrix of data2, 
    sx_tmp = inv(weights{2});
    sx = (dewhiteM{2}*sx_tmp);
%%% calculate the correlation of all componentsts   % match two modality components
    [Corr_matrix]=abs(ica_fuse_corr(mx,sx));
    maxcorr=[];maxcol=[];maxrow=[];
    for j=1:min([ncomps(1),ncomps(2)])
        [mm,ss]=find(Corr_matrix==max(Corr_matrix(:)));
        maxcol(j)=mm(1); maxrow(j)=ss(1); maxcorr(j)=Corr_matrix(mm(1),ss(1));  
        Corr_matrix(mm(1),:)=0;Corr_matrix(:,ss(1))=0;
    end
    [temp, index] = sort(abs(maxcorr));
    index = index(end:-1:1);
  
    if (((step(1)>sparsty_opt_start_step)&&(~stopsign(2)))&&((step(2)>sparsty_opt_start_step)&&(~stopsign(1)))) 
        if ((step(1)>sparsty_opt_start_step)&&(~stopsign(1))&&delta_lambda_hoyer(1)>1e-6)
            if exist('lambda_hoyer1','var')==1    
                [data1,change_flag1(step(1)),sparse_sign1_tmp] = SparsityOpt_HoyerProj(weights{1},data1,1:size(weights{1},1),z_thr,SBR_min_sparstywork,SBR_min(1),step(1),delta_lambda_hoyer(1),sphere_t{1},lambda_hoyer1);
            else
                [data1,change_flag1(step(1)),sparse_sign1_tmp] = SparsityOpt_HoyerProj(weights{1},data1,1:size(weights{1},1),z_thr,SBR_min_sparstywork,SBR_min(1),step(1),delta_lambda_hoyer(1),sphere_t{1});        
            end
            sparse_sign1_clct = [sparse_sign1_clct,sparse_sign1_tmp];
        end
        if ((step(2)>sparsty_opt_start_step)&&(~stopsign(2))&&delta_lambda_hoyer(2)>1e-6)
            if exist('lambda_hoyer2','var')==1    
                [data2,change_flag2(step(2)),sparse_sign2_tmp] = SparsityOpt_HoyerProj(weights{2},data2,1:size(weights{2},1),z_thr,SBR_min_sparstywork,SBR_min(2),step(2),delta_lambda_hoyer(2),sphere_t{2},lambda_hoyer2);
            else
                [data2,change_flag2(step(2)),sparse_sign2_tmp] = SparsityOpt_HoyerProj(weights{2},data2,1:size(weights{2},1),z_thr,SBR_min_sparstywork,SBR_min(2),step(2),delta_lambda_hoyer(2),sphere_t{2});        
            end
            sparse_sign2_clct = [sparse_sign2_clct,sparse_sign2_tmp];
        end
        %%%%reconstruct the data and feed it back for next iteration
       if(step(1)>sparsty_opt_start_step)&&(change_flag1(step(1))>0) 
           icasig1 = data1;
           data1 = inv(weights{1})*data1;
        end

        if(step(2)>sparsty_opt_start_step)&&(change_flag2(step(2))>0) 
           icasig2 = data2;
           data2 = inv(weights{2})*data2;
        end
    end

    % ------------------------- 
    % modifying weights based on correlation between data1 A Matrix and data2 A Matrix  
    weights_beforecorr = weights;  
    if (step(1)>corr_opt_start_step  & step(2)>corr_opt_start_step ) & ( ~stopsign(1) | ~stopsign(2))        
        lossf1_beforecorr(step(1))=Compute_entropy(weights{1},data1_org);
        lossf2_beforecorr(step(2))=Compute_entropy(weights{2},data2_org);

        % A matrix of data1   
        mx = (weights{1}' \ dewhiteM{1}')';
        mx_beforecorr = mx;
        %  A matrix of data2, 
        sx = (weights{2}' \ dewhiteM{2}')';
        sx_beforecorr = sx;
        % calculate the correlation of all componentsts   % match tow modality components
        [Corr_matrix]=abs(ica_fuse_corr(mx,sx));
        maxcorr=[];maxcol=[];maxrow=[];
        for j=1:min([ncomps(1),ncomps(2)])
            [mm,ss]=find(Corr_matrix==max(Corr_matrix(:)));
            maxcol(j)=mm(1); maxrow(j)=ss(1); maxcorr(j)=Corr_matrix(mm(1),ss(1));  
            Corr_matrix(mm(1),:)=0;Corr_matrix(:,ss(1))=0;
        end
        max_corr_sparseopt2_beforecorr = [maxcorr(1)]
        [temp, index] = sort(abs(maxcorr));
        index = index(end:-1:1);
        temp = temp(index);
        ix=find(temp>Connect_threshold);
        if length(ix)>MaxComCon; ix=ix(1:MaxComCon); end
        if ~ isempty(ix)   
            mx_tmp = zeros(size(mx));
            sx_tmp = zeros(size(sx));
            change_flag_sx_corropt1step = 0;
            change_flag_mx_corropt1step = 0;
            for Cons_com=index(1:length(ix)); % constraned componenets                               
                % Updata the weights
                a=mx(:,maxcol(Cons_com))';u=sx(:,maxrow(Cons_com))'; % 
                b=cov(a,u);  b=b(2); tmcorr=b/std(a)/std(u);
                comterm=2*b/var(a)/var(u); coml=length(a);
                
                if ~stopsign(1) %& ~Overindex1                   
                    deltaf=comterm*(u-mean(u)+b*(mean(a)-a)/var(a));% 1st order derivative
                    P=deltaf./norm(deltaf);                         
                    alphak(1)=findsteplength1(-tmcorr^2,-deltaf,a,u,alphak(1),P,0.0001,0.999);
                    aweights_temp=Crate(1)*alphak(1)*P; % 
                    mx_tmp(:,maxcol(Cons_com))=mx(:,maxcol(Cons_com)) + aweights_temp';
                    mx_degree_tmp = acosd(dot(mx_beforecorr(:,maxcol(Cons_com)),mx_tmp(:,maxcol(Cons_com)))/(norm(mx_beforecorr(:,maxcol(Cons_com)),2)*norm(mx_tmp(:,maxcol(Cons_com)),2)));
                    if mx_degree_tmp <= corr_angle_thr
                        mx(:,maxcol(Cons_com)) = mx_tmp(:,maxcol(Cons_com));
                    else
                        %%%%%compute what's the stepsize for correlation direction is when the angle between loading of a component pre-correlation
                        %%%%%and loading of a component after-correlation is larger than the threshold we set(e.g. pi/3)
                        [mx(:,maxcol(Cons_com)),change_flag_mx_tmp] = FindStepSizeAngleExceedThr(mx(:,maxcol(Cons_com)),P',corr_angle_thr);
                        change_flag_mx_corropt1step = change_flag_mx_corropt1step+1;
                        change_flag_mx = [change_flag_mx,change_flag_mx_tmp];                        
                    end
                end
                
                if ~stopsign(2) % & ~Overindex2 
                    deltaf=(comterm*( a-mean(a) + b/var(u)*(mean(u) -u)));% 1st order derivative
                    P=deltaf./norm(deltaf) ;                          
                    alphak(2)=findsteplength1(-tmcorr^2,-deltaf,u,a,alphak(2),P,0.0001,0.999);
                    aweights_temp=Crate(2)*alphak(2)*P; % 
                    sx_tmp(:,maxrow(Cons_com))=sx(:,maxrow(Cons_com)) + aweights_temp'; 
                    sx_degree_tmp = acosd(dot(sx_beforecorr(:,maxrow(Cons_com)),sx(:,maxrow(Cons_com)))/(norm(sx_beforecorr(:,maxrow(Cons_com)),2)*norm(sx(:,maxrow(Cons_com)),2)));
                    if sx_degree_tmp <= corr_angle_thr
                        sx(:,maxrow(Cons_com)) = sx_tmp(:,maxrow(Cons_com));
                    else
                        %%%%%compute what's the stepsize for correlation direction is when the angle between loading of a component pre-correlation
                        %%%%%and loading of a component after-correlation is larger than the threshold we set(e.g. pi/3)
                        [sx(:,maxrow(Cons_com)),change_flag_sx_tmp] = FindStepSizeAngleExceedThr(sx(:,maxrow(Cons_com)),P',corr_angle_thr);
                        change_flag_sx_corropt1step = change_flag_sx_corropt1step+1;
                        change_flag_sx = [change_flag_sx,change_flag_sx_tmp];                       
                    end
                end
            end
            if change_flag_mx_corropt1step>0
                Crate(1) = Crate(1)*0.95;
            end
            if change_flag_sx_corropt1step>0
                Crate(2) = Crate(2)*0.95;
            end
            
            mx_degree_tmp = acosd(dot(mx_beforecorr(:,maxcol(1)),mx(:,maxcol(1)))/(norm(mx_beforecorr(:,maxcol(1)))*norm(mx(:,maxcol(1)))));
            mx_degree = [mx_degree,mx_degree_tmp];
            sx_degree_tmp = acosd(dot(sx_beforecorr(:,maxrow(1)),sx(:,maxrow(1)))/(norm(sx_beforecorr(:,maxrow(1)))*norm(sx(:,maxrow(1)))));
            sx_degree = [sx_degree,sx_degree_tmp];
            
            a=mx(:,maxcol(index(1))); b=sx(:,maxrow(index(1)));
            temp=corrcoef(a,b);  temp=temp(1,2);
            max_corr_corropt_beforeinverse(corropt_step) = temp;%%maximum correlation
            max_corr_monit_beforeinverse = [corropt_step,alphak(1),alphak(2),temp]
            
            if ~stopsign(1)
                temp=weights{1};
                weights{1} = whiteM{1}*mx;   weights{1} =inv(weights{1});
                if max(max(abs(weights{1}))) > MAX_WEIGHT
                    Crate(1)=Crate(1)*0.95;
                    weights{1}=temp;
                end
            end
            
            if  ~stopsign(2)
                temp=weights{2};
                weights{2} = whiteM{2}*sx;
                weights{2} =inv(weights{2});
                if max(max(abs(weights{2}))) > MAX_WEIGHT
                    Crate(2)=Crate(2)*0.95;
                    weights{2}=temp;
                end
            end
            
            sx = (weights{2}' \ dewhiteM{2}')';
            mx = (weights{1}' \ dewhiteM{1}')';
            a=mx(:,maxcol(index(1))); b=sx(:,maxrow(index(1)));
            temp=corrcoef(a,b);  temp=temp(1,2);% correlation
            max_corr_corropt(corropt_step) = temp;%%maximum correlation
            max_corr_monit = [corropt_step,alphak(1),alphak(2),temp]
            if abs(temp)<maxcorr(1);
                disp('Wrong direction !!!! ');    
            end
            lossf3(max(step))=abs(temp);
            oldweights{2} = weights{2};
            oldweights{1} = weights{1};
             corropt_step = corropt_step+1;
        else
            mx_degree_tmp = 0;
            sx_degree_tmp = 0;            
        end      
        
        lossf1_aftercorr(step(1))=Compute_entropy(weights{1},data1_org);
        lossf2_aftercorr(step(2))=Compute_entropy(weights{2},data2_org);
        
    end%--------------------------
    if stopsign(1)==1 & stopsign(2)==1
        laststep=step;
        step=[maxsteps,maxsteps];                % stop when weights stabilize
    end
    
    
end; % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sparse_sign(1) = length(find(sparse_sign1_clct==1));
sparse_sign(2) = length(find(sparse_sign2_clct==1));

if exist('laststep', 'var')
    lrates = lrates(:,1:max(laststep));           % truncate lrate history vector
end

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

function alphak=findsteplength1(fk,deltafk,x1,x2,alphak,P,c1,c2);
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
% end
function [fk_new,change_flag] = FindStepSizeAngleExceedThr(fk,P,angle_thr)
%%%Compute step size for correlation according to preset angle threshold 
a = norm(fk,2);
b = dot(fk,P);
c = norm(P,2);
cos_angle = cos(pi*angle_thr/180);

syms x
eqn = (a*c*cos_angle^2-b^2)*x^2+2*x*(a*b*cos_angle^2-a*b)+(cos_angle^2-1)*a^2== 0;
sol_x = solve(eqn,x);
if isempty(sol_x)
    fk_new = fk;
    change_flag = 0;
else
    if length(sol_x)==2
       sol_x_ind = find(sol_x>0);
       sol_x_f = sol_x(sol_x_ind);
       fk_new = fk+sol_x_f*P;
       angle_check = acosd(dot(fk,fk_new)/(norm(fk,2)*norm(fk_new,2)))
       change_flag = 1;        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Sparsity optimization with Hoyer projection
function [data_tmp,change_flag,sparse_sign,sphere_new] = SparsityOpt_HoyerProj(weights_tmp,data_t,IC_ind,z_thr,SBR_min_sparstywork,SBR_min,step,delta_lambda_hoyer,sphere_tmp,lambda_hoyer)
sparse_sign = 0;
hoyer_small_sign = 0;
s = weights_tmp*data_t;
IC_trimmean = trimmean(s',30);
s_demean = s-repmat(IC_trimmean',1,size(s,2));%%SNP component ususally has non-zero mean,background region may have higher value than signal if mena is skewed, while hoyer projection wouldn't demean first, thus true signal couldn't be enhanced
Hoyer_coeff = HoyerMeasure(s_demean);
change_flag = 0;
data_tmp = data_t;
%%%compute SBR and check if all components achieve the minimum SBR
%%%if not do Hoyer projection to make it achieve the minimum SBR by default
SBR_est = SBR_cal(s_demean,z_thr);
if exist('lambda_hoyer','var')==1 %%if user input a Hoyer index
    if (min(abs(Hoyer_coeff(IC_ind)))< lambda_hoyer)
        sparse_sign = 1;
    elseif (min(abs(Hoyer_coeff(IC_ind))) >= lambda_hoyer)&&(min(SBR_est(IC_ind)) < SBR_min)
        sparse_sign = 1;
        hoyer_small_sign = 1;
    elseif (min(abs(Hoyer_coeff(IC_ind))) >= lambda_hoyer)&&(min(SBR_est(IC_ind)) >= SBR_min)
        sparse_sign = 0;
    end  
else
    if min(SBR_est(IC_ind)) < SBR_min 
        sparse_sign = 1;
    else
        sparse_sign = 0;
    end                
end
%%%identify if user input the hoyer index and the minimum hoyer index of all component is less than the user input hoyer index
Hoyer_coeff_tmp = [];
SBR_tmp = [];
  %%%do hoyer projection
  if step > 5
    if sparse_sign
        change_flag = 1;%%this flag is to monitor the change of the data
        lambda_hoyer_tmp = Hoyer_coeff + delta_lambda_hoyer*ones(size(s_demean,1),1);
        for i=1:length(IC_ind)
            if SBR_est(IC_ind(i))< SBR_min_sparstywork  %%make sure the patterns show up first
                lambda_hoyer_tmp(IC_ind(i))=Hoyer_coeff(IC_ind(i));
            else
                if exist('lambda_hoyer','var')==1
                    if ~hoyer_small_sign %%%%%when users set a larger preset hoyer index, and some components have hoyer index less than the preset hoyer index,then only do hoyer projection on those components
                        if lambda_hoyer_tmp(IC_ind(i)) > lambda_hoyer
                            lambda_hoyer_tmp(IC_ind(i)) = lambda_hoyer;
                        end
                    else %%%%%when users set a preset hoyer index less than the hoyer index of component obtained from data_t, and some components' SBR are less than the threshold,then only do hoyer projection on components with SBR below threshold
                        if SBR_est(IC_ind(i)) > SBR_min 
                            lambda_hoyer_tmp(IC_ind(i))=Hoyer_coeff(IC_ind(i));
                        end                                  
                    end
                else%%%when users don't set a Hoyer index and some components' SBR are less than the threshold,then only do hoyer projection on components with SBR below threshold
                    if SBR_est(IC_ind(i)) > SBR_min
                        lambda_hoyer_tmp(IC_ind(i))=Hoyer_coeff(IC_ind(i));
                    end                                
                end
            end
            
            if (Hoyer_coeff(IC_ind(i))<lambda_hoyer_tmp(IC_ind(i))) 
                Hoyer_coeff_tmp = [Hoyer_coeff_tmp,Hoyer_coeff(IC_ind(i))];
                SBR_tmp = [SBR_tmp,SBR_est(IC_ind(i))];
                s_tmp_l2 = sqrt(sum(s_demean(IC_ind(i),:).^2));%%compute the L2 norm of the current component
                L1s = (sqrt(length(s_demean(IC_ind(i),:)))-(sqrt(length(s_demean(IC_ind(i),:)))-1)*lambda_hoyer_tmp(IC_ind(i)))*s_tmp_l2;                
                s_demean(IC_ind(i),:) = (projfunc(s_demean(IC_ind(i),:)',L1s,s_tmp_l2^2,0))'; %%Given a vector s, find the vector v (not force it to be nonnegative here) having sum(abs(v))=L1s and sum(v.^2)=s_tmp_l2^2 which is closest to s in the euclidian sense.
            end
        end
        data_tmp = s_demean+repmat(IC_trimmean',1,size(s_demean,2));
        if ~isempty(Hoyer_coeff_tmp)
             [delta_lambda_hoyer,sparse_sign,Hoyer_coeff_tmp]
        end
    end
  end
  sphere_new = [];



function Hoyer_coeff = HoyerMeasure(icasig)
%%%%compute Hoyer coefficient for each component
d = size(icasig,2);
Hoyer_coeff = zeros(size(icasig,1),1);
for i = 1:size(icasig,1)
    Hoyer_coeff(i) = (sqrt(d)-norm(icasig(i,:),1)/norm(icasig(i,:),2))/(sqrt(d)-1); 
end

% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SBR_est = SBR_cal(icasig,z_thr)
%%%%compute SBR measure for each component
SBR_est = zeros(size(icasig,1),1);
for i = 1:size(icasig,1)
    cmp_tmp = icasig(i,:);
    cmp_tmp_z = zscore(cmp_tmp);
    [inx]=find((abs(cmp_tmp_z)>=z_thr));
    backgd_vxl_ind = setdiff(1:size(icasig,2),inx);
    SBR_est(i) = mean(abs(cmp_tmp_z(inx)))/(abs(mean(cmp_tmp_z(backgd_vxl_ind)))+std(cmp_tmp_z(backgd_vxl_ind)));
end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Hoyer ptojection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,usediters] = projfunc( s, k1, k2, nn )
% Solves the following problem:
% Given a vector s, find the vector v having sum(abs(v))=k1 
% and sum(v.^2)=k2 which is closest to s in the euclidian sense.
% If the binary flag nn is set, the vector v is additionally
% restricted to being non-negative (v>=0).
%    
% Written 2.7.2004 by Patrik O. Hoyer
%
    
% Problem dimension
N = length(s);

% If non-negativity flag not set, record signs and take abs
if ~nn,
    isneg = s<0;
    s = abs(s);
end

% Start by projecting the point to the sum constraint hyperplane
v = s + (k1-sum(s))/N;

% Initialize zerocoeff (initially, no elements are assumed zero)
zerocoeff = [];

j = 0;
while 1,

    % This does the proposed projection operator
    midpoint = ones(N,1)*k1/(N-length(zerocoeff));
    midpoint(zerocoeff) = 0;
    w = v-midpoint;
    a = sum(w.^2);
    b = 2*w'*v;
    c = sum(v.^2)-k2;
    alphap = (-b+real(sqrt(b^2-4*a*c)))/(2*a);
    v = alphap*w + v;
    
    if all(v>=0),
	% We've found our solution
	usediters = j+1;
	break;
    end
        
    j = j+1;
        
    % Set negs to zero, subtract appropriate amount from rest
    zerocoeff = find(v<=0);
    v(zerocoeff) = 0;
    tempsum = sum(v);
    v = v + (k1-tempsum)/(N-length(zerocoeff));
    v(zerocoeff) = 0;
            
end

% If non-negativity flag not set, return signs to solution
if ~nn,
    v = (-2*isneg + 1).*v;
end

% Check for problems
if max(max(abs(imag(v))))>1e-10,
    error('Somehow got imaginary values!');
end
% end
% end


function [entropy_y] = Compute_entropy(w_,x)
%%%w_ is the weight matrix
%%%x is the data

s_tmp = @(w,x) w*x;
w0 = @(w,x) zeros(size(s_tmp(w,x)));
u = @(w,x) s_tmp(w,x)+w0(w,x);
logy = @(w,x) -log((1+exp(-u(w,x))));% + repmat(1000*eps*sign(rand(size(u(w,x),1),1)-.5),1,size(u(w,x),2));
logDv = @(w,x) logy(w,x) + log(1-exp(logy(w,x)));
Wn = @(w) w*w';
lambda = @(w) sqrt(eig(Wn(w)));

remove = isinf(sum(abs(logDv(w_,x)))) | isnan(sum(abs(logDv(w_,x))));
x = x(:,~remove);%%%%discard those infinity or NaN voxels from x

ica.func = @(w,x) (mean(sum(abs(logDv(w,x)))) + sum(log(abs(lambda(w)))));

% V = size(x,2);
% s_tmp = @(w,x) w*x;
% w0 = @(w,x) ones(size(s_tmp(w,x)));
% u = @(w,x) s_tmp(w,x)+w0(w,x);
% logy = @(w,x) -log((1+exp(-u(w,x)))) + repmat(1000*eps*sign(rand(size(u(w,x),1),1)-.5),1,size(u(w,x),2));%%% +repmat(1000*eps*sign(rand(size(u(w,x),1),1)-.5),1,size(u(w,x),2)) to avoid exp(logy(w,x)) to be infinity
% logDv = @(w,x) logy(w,x) + log(1-exp(logy(w,x)));
% Wn = @(w) w*w';
% lambda = @(w) sqrt(eig(Wn(w)));
% ica.func = @(w,x) sum(sum(abs(logDv(w,x)))) + V*sum(log(abs(lambda(w))));

entropy_y = ica.func(w_,x);