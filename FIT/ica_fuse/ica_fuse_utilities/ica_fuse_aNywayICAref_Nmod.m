% aNywayICAref_Nmod() - Perform reference-constrained aNyway independent component analysis (aNyway ICA with reference) decomposition,
% aNyway ICA with reference maximize the independence of components for each modality based on infomax ICA algorithm (Bell & Sejnowski (1995)) 
% and miminize the mutual information among subspace component vectors(i.e., maximize the cross-correlation of loadings among modalities and reference) 
% based on Gaussian independent vector analysis (IVA-G, Anderson, Adali & Li (2012)) via a shared weight matrix model without orthogonality constraints
% Usage:[weights,sphere,activations,A] = aNywayICAref_Nmod(data,'ncomps',ncomps_num,'iva_lambda',lambda_ivaRef,'reference',reference);
% Input_Variable:
%
% data        = cell type input data:{data1, data2,... dataN}; for each dataset, the dimension can be subject-by-features (sMRI) or subject*timepoints-by-features (fMRI).
%               Note: If data consists of multiple discontinuous epochs,
%               each epoch should be separately baseline-zero'd using:
%                  >> data = rmbase(data,frames,basevector);
% Key Keywods value
% 'ncomps' = an array contains values of components for each modality (aNyway ICA allows different component numbers for different modalities)
% 'lambda_ivaRef' = a value between 0 and 1, the regularizer to balance between independence and correlation
% 'reference' = a vector (subject-by-1) as the reference constraint (e.g., symptom score, cognitive performance, polygenic risk score)

% Optional_Keywords       Keyword_Values                  Default_Values
% 'pca'       = [N] decompose a principal component     (default -> 0=off)
%               subspace of the data. Value is the number of PCs to retain.
% 'sphering'  = ['on'/'off'] flag sphering of data      (default -> 'off')
% 'weights'   = [W] initial weight matrix               (default -> eye())
%                            (Note: if 'sphering' 'off', default -> spher())
% 'lrate'     = [rate] initial ICA learning rate (<< 1) (default -> heuristic)
% 'block'     = [N] ICA block size (<< datalength)      (default -> heuristic)
% 'anneal'    = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended controls speed of convergence)                     
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
% 'iva_lambda_prealign' = a value between 0 and 1 to change the weight update factor for IVA ref prealignment (to avoid weight blowup)
%
% Output_Variables [RO = output in reverse order of projected mean variance
%                        unless starting weight matrix passed ('weights' above)]
% weights     = ICA weight matrix (comps,chans)     [RO]
% sphere      = data sphering matrix (chans,chans) = spher(data)
%               Note: unmixing_matrix = weights*sphere {sphering off -> eye(chans)}
% activations = activation time courses of the output components (ncomps,frames*epochs)
% A           = loading matrix (subject,ncomps)

%% Toolbox Citation:
% Duan, K., Silva, R. F., Liu, J., Agcaoglu, O., & Calhoun, V. D. (2023, April). Any-Way Independent Component Analysis with Reference. In 2023 IEEE 20th International Symposium on Biomedical Imaging (ISBI) (pp. 1-4). IEEE.
% Duan, K., Silva, R. F., Liu, J., & Calhoun, V. D. (2020, July). aNy-way independent component analysis. In 2020 42nd Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC) (pp. 1770-1774). IEEE.

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


function [weights,sphere,activations,A_f,bias,signs,lrates,y] = ica_fuse_aNywayICAref_Nmod(data,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)

if nargin < 2
    help aNywayICAref_Nmod
    return
end

%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      % anneal by multiplying lrate by this original 0,9 to 0.95 changed by JL
DEFAULT_EXTANNEAL    = 0.98;      % or this if extended-ICA
DEFAULT_MAXSTEPS     = 1000;      %%DEFAULT_MAXSTEPS     = 512;     % stop training after this many steps 
DEFAULT_MOMENTUM     = 0.0;       % default momentum weight

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       %%% when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.95;       % %%DEFAULT_RESTART_FAC  = 0.9;  if weights blowup, restart with lrate
% lower by this factor
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
% % DEFAULT_BLOCK_NEW = DEFAULT_BLOCK;

% - may need adjustment!
% Extended-ICA option:
DEFAULT_EXTENDED     = 0;         % default off
DEFAULT_EXTBLOCKS    = 1;         % number of blocks per kurtosis calculation
DEFAULT_NSUB         = 1;         % initial default number of assumed sub-Gaussians
% for extended-ICA
DEFAULT_EXTMOMENTUM  = 0.5;       % momentum term for computing extended-ICA kurtosis
MAX_KURTSIZE         = 6000;      % max points to use in kurtosis calculation
MIN_KURTSIZE         = 20001;      % minimum good kurtosis size (flag warning)
SIGNCOUNT_THRESHOLD  = 25;        % raise extblocks when sign vector unchanged
% after this many steps
SIGNCOUNT_STEP       = 2;         % extblocks increment factor

DEFAULT_SPHEREFLAG   = 'none';    %%do not use it for aNy-way ICA %%DEFAULT_SPHEREFLAG   = 'on';  % use the sphere matrix as the default
%   starting weight matrix
DEFAULT_PCAFLAG      = 'off';     % don't use PCA reduction
DEFAULT_POSACTFLAG   = 'off';     %%do not use it for aNy-way ICA %%DEFAULT_POSACTFLAG   = 'on';% use posact()
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule

DEFAULT_IVA_maxstep = 50; %%%run IVA with reference 50 steps for alignment
DEFAULT_LAMBDA_IVAref = 0.15;
IVA_weight_adjust_start_step = 80;
DEFAULT_LAMBDA_IVAref_ANEAL_FAC = 0.95;
DEFAULT_reference = ones(size(data{1},1),1);
DEFAULT_IVAREFPREALIGN_LAMBDA_ADJUST_FAC = 1/150;%%%default adjust factor for IVA with reference

lambda_IVAref = DEFAULT_LAMBDA_IVAref;
lambda_IVAref_aneal_fac = DEFAULT_LAMBDA_IVAref_ANEAL_FAC;
ivaRefprealign_lambda_adjust_fac = DEFAULT_IVAREFPREALIGN_LAMBDA_ADJUST_FAC;

ENDURANCE = -1e-3; %%the maximumlly allowed descending trend of entropy;

%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
if (nargin> 1 & rem(nargin,2) == 0)
    fprintf('runica() : Even number of input arguments???')
    return
end
for i = 3:2:nargin % for each Keyword
    Keyword = eval(['p',int2str((i-3)/2 +1)]);
    Value = eval(['v',int2str((i-3)/2 +1)]);
    if ~isstr(Keyword)
        fprintf('runica() : keywords must be strings')
        return
    end
    Keyword = lower(Keyword); % convert upper or mixed case to lower
    
    if strcmp(Keyword,'weights') | strcmp(Keyword,'weight')
        if isstr(Value)
            fprintf(...
                'runica() : weights value must be a weight matrix or sphere')
            return
        else
            weights = Value;
            wts_passed =1;
        end
    elseif strcmp(Keyword,'ncomps')
        if isstr(Value)
            fprintf('runica() : ncomps value must be an integer')
            return
        else
            %%%added for aNy-way ICA with reference
            ncomps = Value;
            urchans = ncomps;
        end
        if ncomps < urchans & ncomps ~= Value
            fprintf('runica() : Use either PCA or ICA dimension reduction');
            return
        end  
        ncomps = Value;
        if ~ncomps,
            ncomps = chans;
        end
    elseif strcmp(Keyword,'pca')
        if ncomps < urchans & ncomps ~= Value
            fprintf('runica() : Use either PCA or ICA dimension reduction');
            return
        end
        if isstr(Value)
            fprintf(...
                'runica() : pca value should be the number of principal components to retain')
            return
        end
        pcaflag = 'on';
        ncomps = Value;
        if ncomps >= chans | ncomps < 1,
            fprintf('runica() : pca value must be in range [1,%d]\n',chans-1)
            return
        end
        chans = ncomps;
    elseif strcmp(Keyword,'posact')
        if ~isstr(Value)
            fprintf('runica() : posact value must be on or off')
            return
        else
            Value = lower(Value);
            if ~strcmp(Value,'on') & ~strcmp(Value,'off'),
                fprintf('runica() : posact value must be on or off')
                return
            end
            posactflag = Value;
        end
    elseif strcmp(Keyword,'lrate')
        if isstr(Value)
            fprintf('runica() : lrate value must be a number')
            return
        end
        lrate = Value;
        if lrate>MAX_LRATE | lrate <0,
            fprintf('runica() : lrate value is out of bounds');
            return
        end
        if ~lrate,
            lrate = DEFAULT_LRATE;
        end
    elseif strcmp(Keyword,'iva_lambda')
        if isstr(Value)
            fprintf('runica() : iva_lambda value must be a number')
            return
        end
        lambda_IVAref = Value;
        if lambda_IVAref>1 | lambda_IVAref <0,
            fprintf('runica() : lambda_IVAref is out of bounds');
            return
        end
        if ~lambda_IVAref,
            lambda_IVAref = DEFAULT_LAMBDA_IVAref;
        end
        
    elseif strcmp(Keyword,'iva_lambda_prealign') 
        if isstr(Value)
            fprintf('runica() : iva_lambda_prealign value must be a number')
            return
        end
        ivaRefprealign_lambda_adjust_fac = Value;
        if ivaRefprealign_lambda_adjust_fac>1 | ivaRefprealign_lambda_adjust_fac <0,
            fprintf('runica() : iva_lambda_prealign value is out of bounds');
            return
        end
        if ~ivaRefprealign_lambda_adjust_fac,
            ivaRefprealign_lambda_adjust_fac = DEFAULT_IVAREFPREALIGN_LAMBDA_ADJUST_FAC;
        end               
    elseif strcmp(Keyword,'reference')
        if isstr(Value)
            fprintf('runica(): reference value must be a vector')
            return
        elseif(length(Value)~=size(data{1},1))
            fprintf('runica(): length of reference should be subject number')
            return        
        end
        reference = Value;
        if ~reference,
            reference = DEFAULT_reference;
        end    
    elseif strcmp(Keyword,'ivamaxstep')
        if isstr(Value)
            fprintf('runica(): IVA_maxstep must be a integer')
            return      
        end
        IVA_maxstep = Value;
        if ~IVA_maxstep,
            IVA_maxstep = DEFAULT_IVA_maxstep;
        end

    elseif strcmp(Keyword,'block') | strcmp(Keyword,'blocksize')
        if isstr(Value)
            fprintf('runica() : block size value must be a number')
            return
        end
        block = Value;
        if ~block,
% %             block = DEFAULT_BLOCK;
            block = DEFAULT_BLOCK_NEW;
        end
    elseif strcmp(Keyword,'stop') | strcmp(Keyword,'nochange') ...
            | strcmp(Keyword,'stopping')
        if isstr(Value)
            fprintf('runica() : stop wchange value must be a number')
            return
        end
        nochange = Value;
    elseif strcmp(Keyword,'maxsteps') | strcmp(Keyword,'steps')
        if isstr(Value)
            fprintf('runica() : maxsteps value must be an integer')
            return
        end
        maxsteps = Value;
        if ~maxsteps,
            maxsteps   = DEFAULT_MAXSTEPS;
        end
        if maxsteps < 0
            fprintf('runica() : maxsteps value must be a positive integer')
            return
        end
    elseif strcmp(Keyword,'anneal') | strcmp(Keyword,'annealstep')
        if isstr(Value)
            fprintf('runica() : anneal step constant must be a number (0,1)')
            return
        end
        annealstep = Value;
        if annealstep <=0 | annealstep > 1,
            fprintf('runica() : anneal step value must be (0,1]')
            return
        end
    elseif strcmp(Keyword,'annealdeg') | strcmp(Keyword,'degrees')
        if isstr(Value)
            fprintf('runica() : annealdeg value must be a number')
            return
        end
        annealdeg = Value;
        if ~annealdeg,
            annealdeg = DEFAULT_ANNEALDEG;
        elseif annealdeg > 180 | annealdeg < 0
            fprintf('runica() : annealdeg value is out of bounds [0,180]')
            return            
        end
    elseif strcmp(Keyword,'momentum')
        if isstr(Value)
            fprintf('runica() : momentum value must be a number')
            return
        end
        momentum = Value;
        if momentum > 1.0 | momentum < 0
            fprintf('runica() : momentum value is out of bounds [0,1]')
            return
        end
    elseif strcmp(Keyword,'sphering') | strcmp(Keyword,'sphereing') ...
            | strcmp(Keyword,'sphere')
        if ~isstr(Value)
            fprintf('runica() : sphering value must be on, off, or none')
            return
        else
            Value = lower(Value);
            if ~strcmp(Value,'on') & ~strcmp(Value,'off') & ~strcmp(Value,'none'),
                fprintf('runica() : sphering value must be on or off')
                return
            end
            sphering = Value;
        end
    elseif strcmp(Keyword,'bias')
        if ~isstr(Value)
            fprintf('runica() : bias value must be on or off')
            return
        else
            Value = lower(Value);
            if strcmp(Value,'on')
                biasflag = 1;
            elseif strcmp(Value,'off'),
                biasflag = 0;
            else
                fprintf('runica() : bias value must be on or off')
                return
            end
        end
    elseif strcmp(Keyword,'specgram') | strcmp(Keyword,'spec')
        if isstr(Value)
            fprintf('runica() : specgram argument must be a vector')
            return
        end
        srate = Value(1);
        if (srate < 0)
            fprintf('runica() : specgram srate must be >=0')
            return
        end
        if length(Value)>1
            loHz = Value(2);
            if (loHz < 0 | loHz > srate/2)
                fprintf('runica() : specgram loHz must be >=0 and <= srate/2')
                return
            end
        else
            loHz = 0; % default
        end
        if length(Value)>2
            hiHz = Value(3);
            if (hiHz < loHz | hiHz > srate/2)
                fprintf('runica() : specgram hiHz must be >=loHz and <= srate/2')
                return
            end
        else
            hiHz = srate/2; % default
        end
        if length(Value)>3
            Hzinc = Value(4);
            if (Hzinc<=0 | Hzinc>hiHz-loHz)
                fprintf('runica() : specgram Hzinc must be >0 and <= hiHz-loHz')
                return
            end
        else
            Hzinc = 1; % default
        end
        if length(Value)>4
            Hzframes = Value(5);
            if (Hzframes<0 | Hzframes > size(data1,2))
                fprintf('runica() : specgram frames must be >=0 and <= data length')
                return
            end
        else
            Hzframes = size(data1,2); % default
        end
    elseif strcmp(Keyword,'extended') | strcmp(Keyword,'extend')
        if isstr(Value)
            fprintf('runica() : extended value must be an integer (+/-)')
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
                        ' runica() warning: kurtosis values inexact for << %d points.\n',...
                        MIN_KURTSIZE);
                end
            end
        end
    elseif strcmp(Keyword,'verbose')
        if ~isstr(Value)
            fprintf('runica() : verbose flag value must be on or off')
            return
        elseif strcmp(Value,'on'),
            verbose = 1;
        elseif strcmp(Value,'off'),
            verbose = 0;
        else
            fprintf('runica() : verbose flag value must be on or off')
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
        fprintf('runica() : unknown flag')
        return
    end
end
%

% %-----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%Data initialization and PCA
data_org = data; %%save the original data for later usage.
subj_N = size(data{1},1);
data_totalvar = zeros(length(data),1);
Var_retained_data_org = zeros(length(data),1);

%%%%PCA on each modality to reduce the dimension 
for m = 1:length(data)
    c = ncomps(m);
    data_demean_row{m} = data{m}-mean(data{m},2);%%remove row mean.
    data_demean_row_column{m} = data_demean_row{m}-mean(data_demean_row{m},1);%%then remove column mean  
    %%%%%compute the covariance of principal components
    if size(data_demean_row_column{m},1) <= size(data_demean_row_column{m},2) %%%short matrices
        data_cov{m} = data_demean_row_column{m}*data_demean_row_column{m}'/(size(data_demean_row_column{m},2)-1);%%sample variance
        [U1_data{m},S1_data{m}] = eigs(data_cov{m},ncomps(m)); %%%Do SVD on demeaned data for dimensiom reduction,eigen value is in ascending order
        data_totalvar(m) = trace(data_cov{m}); 
        V1_transp{m} = pinv(sqrt(S1_data{m}))*(U1_data{m}(:,1:c)'*data_demean_row_column{m});%%Cm-by-v dimension
        V1{m} = V1_transp{m}';
    else   %%%tall matrices
        data_cov{m} = data_demean_row_column{m}'*data_demean_row_column{m}/(size(data_demean_row_column{m},1)-1);
        [V1{m},S1_data{m}] = eigs(data_cov{m},ncomps(m)); 
        data_totalvar(m) = trace(data_cov{m}); 
        U1_data{m} = (data_demean_row_column{m}*V1{m})*pinv(sqrt(S1_data{m}));
        V1_transp{m} = V1{m}';
    end

    KInv{m} = std(V1{m});
    K{m} = 1./KInv{m};   
    KtimesV1_transp{m} = repmat(K{m}',1,size(V1_transp{m},2)).*V1_transp{m};
    Var_retained_data_org(m) = 100*(trace(S1_data{m})/data_totalvar(m));
    fprintf('For %d-th modality, final component number is %d, and %g%% of (non-zero) variance retained.\n',m,c,Var_retained_data_org(m));

    whiteM{m} = repmat(K{m}',1,size(U1_data{m},1)).*(U1_data{m}(:,1:c))';
    dewhiteM{m} = pinv(whiteM{m});
    V1_transp_inv{m} = pinv(V1_transp{m});   
    VdevideByK{m} = V1_transp_inv{m}.*repmat(KInv{m},size(V1{m},1),1);   
    data{m} = KtimesV1_transp{m};%%%%data after PCA.
        
    [chans(m) frames(m)] = size(data{m}); % determine the data size
    urchans(m) = chans(m);  % remember original data channels
    datalength(m) = frames(m);
end

%%%re-order the principal component according to its correlation with the reference to provide a good starting point for IVA with reference 
[scv_num,ICnum_min_mod_ind] = min(ncomps);
mod_moreIC_ind = setdiff(1:length(ncomps),ICnum_min_mod_ind);
for i = 1:length(data)
    mod_refcorr_abs{i} = abs(corr(reference,dewhiteM{i}));
    [mod_refcorr_abs_max{i},mod_refcorr_abs_max_ind{i}] = sort(mod_refcorr_abs{i},'descend'); 

    %%%%%reorder the principal components
    V1_transp{i} = V1_transp{i}(mod_refcorr_abs_max_ind{i},:); 
    V1{i} = V1_transp{i}';
    KInv{i} = std(V1{i});
    K{i} = 1./KInv{i};   
    KtimesV1_transp{i} = repmat(K{i}',1,size(V1_transp{i},2)).*V1_transp{i};
    V1_transp_inv{i} = pinv(V1_transp{i});   
    whiteM{i} = whiteM{i}(mod_refcorr_abs_max_ind{i},:);
    dewhiteM{i} = pinv(whiteM{i});
    VdevideByK{i} = V1_transp_inv{i}.*repmat(KInv{i},size(V1{i},1),1);
    data{i} = KtimesV1_transp{i};
end 

%-----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
DEFAULT_LRATE        = max(0.01*(0.015./log(chans)))*ones(1,length(chans));%%% DEFAULT_LRATE        = 0.015./log(chans);
% heuristic default - may need adjustment
%   for large or tiny data sets!
DEFAULT_BLOCK        = floor(sqrt(frames/3));  % heuristic default %% a vector stores the block size for each modality
block_scale = 10;
DEFAULT_BLOCK_NEW = block_scale.*DEFAULT_BLOCK;
%%%May 4, 2021, add labmda for refrence
lambda_ref = 1;
%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargout < 2,
    fprintf('runica()  - needs at least two output arguments.\n');
    return
end
epochs = 1;							 % do not care how many epochs in data

pcaflag    = DEFAULT_PCAFLAG;
sphering   = DEFAULT_SPHEREFLAG;     % default flags
posactflag = DEFAULT_POSACTFLAG;
verbose    = DEFAULT_VERBOSE;
% block      = DEFAULT_BLOCK;          % heuristic default - may need adjustment!
block      = DEFAULT_BLOCK_NEW;          % heuristic default - may need adjustment!
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
kurtsize   = MAX_KURTSIZE*ones(1,length(data));
signsbias  = 0.02;                   % bias towards super-Gaussian components
extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
nsub       = DEFAULT_NSUB;
wts_blowup = zeros(1,length(data));  % flag =1 when weights too large
wts_passed = 0;                      % flag weights passed as argument

trendPara  = ENDURANCE; %depends on the requirement on connection; the more negative,the stronger the contrains ,that may cause overfitting

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
    fprintf(' runica(): number of components must be 1 to %d.\n',chans);
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
                ' runica(): weight matrix must have %d rows, %d columns.\n', ...
                chans,ncomps);
            return;
        end
    end
end;
%
%%%%%%%%%%%%%%%%%%%%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if frames<chans,
    fprintf(' runica(): data length %d < data channels %f!\n',frames,chans)
    return
elseif block < 2,
    fprintf(' runica(): block size %d too small!\n',block)
    return
elseif block > frames,
    fprintf(' runica(): block size exceeds data length!\n');
    return
elseif floor(epochs) ~= epochs,
    fprintf(' runica(): data length is not a multiple of the epoch length!\n');
    return
elseif nsub > ncomps
    fprintf(' runica(): there can be at most %d sub-Gaussian components!\n',ncomps);
    return
end;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
    for i = 1:length(data)
            fprintf( ...
        '\nInput data size [%d,%d] = %d channels, %d frames.\n', ...
        chans(i),frames(i),chans(i),frames(i));
    end
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
    for i = 1:length(data)
       fprintf('Data %d, Initial learning rate will be %g, block size %d.\n',i, lrate(i),block(i)); 
    end
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
    for i = 1:length(data)
        fprintf('Final training data{1} range: %g to %g\n', ...
        min(min(data{i})),max(max(data{i})));
    end  
end

%
%%%%%%%%%%%%%%%%%%% Perform PCA reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Reducing the data to %d principal dimensions...\n',ncomps);
    for i = 1:length(data)
        [eigenvectors{i},eigenvalues{i},data_whitened{i}] = pcsquash(data_org{i},ncomps(i));
    end
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
    sphere = [];
    for i = 1:length(data)
        sphere{i} = 2.0*inv(sqrtm(cov(data{i}'))); % find the "sphering" matrix = spher()
    end

    if ~weights,
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
        end
        weights=[];
        for i = 1:length(data)
           weights{i} = eye(ncomps(i),chans(i)); % begin with the identity matrix
        end
    else % weights given on commandline
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
        end
    end
    if verbose,
        fprintf('Sphering the data ...\n');
    end
    
    for i = 1:length(data)
       data{i} = sphere{i}*data{i};      % actually decorrelate the electrode signals
    end
elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~weights
        if verbose,
            fprintf('Using the sphering matrix as the starting weight matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        weights = [];
        sphere = [];
        for i = 1:length(data)
            sphere{i} = 2.0*inv(sqrtm(cov(data{i}'))); % find the "sphering" matrix = spher()
            weights{i} = eye(ncomps(i),chans(i))*sphere{i}; % begin with the identity matrix
            sphere{i} = eye(chans(i));                 % return the identity matrix  
        end      
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        for i = 1:length(data)
            sphere{i} = eye(chans(i));                 % return the identity matrix  
        end
    end
elseif strcmp(sphering,'none')
    for i = 1:length(data)
        sphere{i} = eye(chans(i));                 % return the identity matrix  
    end
    if ~weights
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        weights=[];      
        for i = 1:length(data)
            weights{i} = eye(ncomps(i),chans(i));                 % begin with the identity matrix
        end
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
    end
    for i = 1:length(data)
        sphere{i} = eye(chans(i));                 % return the identity matrix  
    end
    if verbose,
        fprintf('Returned variable "sphere" will be the identity matrix.\n');
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
lastt=fix((datalength./block-1).*block+1);
changes = [];
degconst = 180./pi;
startweights = weights;
prevweights = startweights;
oldweights = startweights;
lrates = zeros(length(data),maxsteps);
delta_concat = zeros(1,sum(chans.*ncomps));

BI = []; delta = []; prevwtchange = []; oldwtchange = [];
onesrow = []; bias = []; signs = [];
for i = 1:length(data)
    BI{i}=block(i)*eye(ncomps(i),ncomps(i));
    delta{i}=zeros(1,chans(i)*ncomps(i));
    prevwtchange{i} = zeros(chans(i),ncomps(i));
    oldwtchange{i} = zeros(chans(i),ncomps(i));
    onesrow{i} = ones(1,block(i));
    bias{i} = zeros(ncomps(i),1);
    signs{i} = ones(1,ncomps(i));    % initialize signs to nsub -1, rest +1
    for k=1:nsub
        signs{i}(k) = -1;
    end
end

if extended & extblocks < 0 & verbose,
    fprintf('Fixed extended-ICA sign assignments:  ');
    for i = 1:length(data)
        for k=1:ncomps
            fprintf('%d ',signs{i}(k));
        end; fprintf('\n');
    end
end

for i = 1:length(data)
    signs{i} = diag(signs{i}); % make a diagonal matrix
    oldsigns{i} = zeros(size(signs{i}));
    old_kk{i} = zeros(1,ncomps(i));   % for kurtosis momemtum
end

change = zeros(1,length(data));
signcount = zeros(1,length(data));     % counter for same-signs
signcounts = [];
urextblocks = extblocks;    % original value, for resets
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
step = zeros(1,length(data));
blockno = ones(1,length(data));  % running block counter for kurtosis interrupts
stopsign = zeros(1,length(data));      
angledelta = zeros(1,length(data));
index1 = zeros(1,length(data));
                
p = 1:min(ncomps);
K_scv = min(ncomps); 
batch_ind_max_tmp = [];
for i = 1:length(data)
    ICA_cost{i} = zeros(maxsteps,1);
    ICA_cost2{i} = zeros(maxsteps,1);
    batch_ind{i} = 1:block(i):lastt(i);
    batch_ind_max_tmp = [batch_ind_max_tmp,max(batch_ind{i})];
end
ICAcost_threemod = zeros(maxsteps,1);
IVA_cost = zeros(maxsteps,1);
batch_ind_max = max(batch_ind_max_tmp);

batch_num = zeros(length(data),1);
for i = 1:length(data)
    batch_num(i) = length(batch_ind{i});
end
[batch_num_max,batch_num_ind] = max(batch_num(:));
batch_num_missed = batch_num_max-batch_num(:);
smaller_batch_num_modlty_ind = find(batch_num_missed>0);
max_batch_num_modlty_ind = setdiff(1:length(data),smaller_batch_num_modlty_ind);
batch_ind_new = cell(length(data),1);
for j = 1:length(max_batch_num_modlty_ind)
    batch_ind_new{max_batch_num_modlty_ind(j)} = batch_ind{max_batch_num_modlty_ind(j)};
end
data_pca = data;
%-------------------------------------------------------------------------------------------------
%%%%%%%Apply IVA with reference constraint to pre-align the SCVs to provide a good starting point for aNy-way ICA with reference
fprintf('Beginning IVA with reference training for prealignment.\n');
fprintf('IVA with reference prealignment lambda is %f.\n',ivaRefprealign_lambda_adjust_fac);
fprintf('aNy-way ICA with reference, IVA with reference regularizer lambda is %f.\n',lambda_IVAref);
[weights_ivaRef_new,iva_stepF,IVA_cost] = ica_fuse_IVAref_Nway(data,data_org,reference,VdevideByK,ivaRefprealign_lambda_adjust_fac,dewhiteM);

%%%%Project the data to the weight matrix obtained from IVA with reference
for m = 1:length(data)
    data{m} = weights_ivaRef_new{m}*data{m};
    VdevideByK{m} = VdevideByK{m}/weights_ivaRef_new{m};
    whiteM{m} = weights_ivaRef_new{m}*whiteM{m};
    dewhiteM{m} = dewhiteM{m}/weights_ivaRef_new{m};
end
%-------------------------------------------------------------------------------------------------

%%-------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%aNy-way ICA with reference constraint%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Beginning aNyway ICA with reference training ...\n');
weights_ifmx = cell(length(data),1);
entropy = zeros(length(data),1);
entropychange = zeros(length(data),1);
while ~isempty(find(step < maxsteps))
    for d = 1:length(data)
        ICA_cost{d}(step(d)+1) = Compute_entropy_infmx(weights{d},data{d},bias{d});
        entropy(d) = ICA_cost{d}(step(d)+1);
    end
    if ~isempty(find(stopsign==0))  
        %%%%%%%%%%%%%optimize aNy-way ICA with reference using mini batch              
        if ~isempty(find(step>=0))
            if unique(step)>150
                stopsign
            end
            IVA_weights_factor = 0.05*sum(block);
            for h = 1:length(data)
                permuteVec{h} = randperm(datalength(h)); %shuffle features at each step
            end
            
            %%%expand those modalities with smaller batch number to make all modalities have same batch numbers
            if ~ isempty(smaller_batch_num_modlty_ind)           
                for mis = 1:length(smaller_batch_num_modlty_ind)
                    batch_duplct = [];
                    while (length(batch_ind{smaller_batch_num_modlty_ind(mis)})+length(batch_duplct))< batch_num_max
                        smaller_batch_modlty_perm = randperm(batch_num(smaller_batch_num_modlty_ind(mis)));
                        if batch_num_missed(smaller_batch_num_modlty_ind(mis)) <= length(batch_ind{smaller_batch_num_modlty_ind(mis)})
                            batch_duplct = [batch_duplct, batch_ind{smaller_batch_num_modlty_ind(mis)}(smaller_batch_modlty_perm(1:batch_num_missed(smaller_batch_num_modlty_ind(mis))))];                     
                        else
                            batch_duplct = [batch_duplct, batch_ind{smaller_batch_num_modlty_ind(mis)}(smaller_batch_modlty_perm(1:length(batch_ind{smaller_batch_num_modlty_ind(mis)})))];
                            batch_num_missed(smaller_batch_num_modlty_ind(mis))= batch_num_missed(smaller_batch_num_modlty_ind(mis))-length(batch_ind{smaller_batch_num_modlty_ind(mis)});
                        end
                    end
                    batch_ind_new{smaller_batch_num_modlty_ind(mis)} = [batch_ind{smaller_batch_num_modlty_ind(mis)},batch_duplct];
                end            
            end
            IVAjoint_func = 0;
            IVAMarginalEntropy_func = 0;
            ICAcost_threemod_sum = 0;
            for b = 1:batch_num_max
                %%%%1,for each modality,compute infomax weight update
                ICAcost_threemod = zeros(length(data),1);
                for i = 1:length(data)            
                    IVAjoint_diff_weights{i} = zeros(K_scv,K_scv);
                    IVAmarginal_diff_weights{i} = zeros(K_scv,K_scv);
                    IVA_diff_weights{i} = eye(size(weights{i}));
                    S_inv_deriv_seprt{i} = zeros(frames(i),K_scv);
                    weights_ifmx{i} = zeros(size(weights{i}));
                    weights_ICA{i} = zeros(size(weights{i}));
                    if (~stopsign(i))
                        SInv{i} = VdevideByK{i}/weights{i};%%%Update inverse of S with updated weight matrix                    
                        A{i} = data_org{i}*SInv{i};
                        if biasflag            
                            u=weights{i}*data{i}(:, permuteVec{i}(batch_ind_new{i}(b):batch_ind_new{i}(b)+block(i)-1)) + bias{i}*onesrow{i};
                        else
                            u=weights{i}*data{i}(:, permuteVec{i}(batch_ind_new{i}(b):batch_ind_new{i}(b)+block(i)-1));
                        end
                        if ~extended
                            %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%                
                            y=1./(1+exp(-u));
                            weights_ifmx{i} = (BI{i}+(1-2*y)*u')*weights{i};
      %                     weights{1}=weights{1}+lrate(1)*(BI{1}+(1-2*y)*u')*weights{1};                   %                
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else % extended-ICA
                            %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
                            y=tanh(u); 
                            weights_ifmx{i} = (BI{i}-signs{i}*y*u'-u*u')*weights{i};%
                            %weights{1} = weights{1} + lrate(1)*(BI{1}-signs{1}*y*u'-u*u')*weights{1};          %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
                        if biasflag
                            if ~extended
                                %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%                    
                                bias{i} = bias{i} + lrate(i)*sum((1-2*y)')'; % for logistic nonlin. %
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            else % extended
                                %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                bias{i} = bias{i} + lrate(i)*sum((-2*y)')';  % for tanh() nonlin.   %
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            end
                        end

                        %%%check whether ICA weight update yields to blowup or not, if yes, do not do IVA ref update, and reset the initial weight matrix and reduce the learning rate to restart again
                        weights_ICA{i} = weights{i}+lrate(i)*weights_ifmx{i};                      
                        ICAcost_threemod(i) = Compute_entropy_infmx(weights_ICA{i},data{i},bias{i});
                        if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             weights_ICA{i} = weights_ICA{i} + momentum*prevwtchange{i};                
                             prevwtchange{i} = weights_ICA{i}-prevweights{i};                      
                             prevweights{i} = weights_ICA{i};                                  
                        end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if max(max(abs(weights_ICA{i}))) > MAX_WEIGHT
                             wts_blowup(i) = 1;
                        end

                        %%%%2,for each modality,compute IVA ref weight update
                        %%%------------------------2.1,IVA ref joint entropy weight update------------------------------        
                        %%%%%%%%%%%%%%%%%%%%IVA weights update starts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if (~wts_blowup(i)) 
                            B{i} = SInv{i}(:,p);
                            [V2{i},S2{i}] = eig(B{i}'*B{i});
                            S2{i} = sqrt(diag(S2{i}));
                            U2{i} = B{i}/((V2{i})'.*repmat(S2{i},1,size(V2{i},2)));
                            IVAjoint_func = IVAjoint_func-sum(log(abs(S2{i}(:))));
                            Q{i} = zeros(K_scv,K_scv);
                            Q{i}(p,:) = V2{i};

                            %-------------------------------------------
                            F{i} = U2{i}'*VdevideByK{i}(:,p); 
                            FdevideW{i} = F{i}/oldweights{i}(p,p);
                            T{i} = FdevideW{i}*Q{i};
                            QdevideT{i} = Q{i}/T{i};   
                            %%%%-----------natural gradient(regular gradient*W^T*W)
                            IVAjoint_diff_weights{i} = (oldweights{i}(p,p)\((QdevideT{i})*F{i}))'*oldweights{i}(p,p);  
                        end
                        %%------------------------2.1,IVA joint entropy weight update------------------------------                    
                     end             
                end
                ICAcost_threemod_sum = ICAcost_threemod_sum+sum(ICAcost_threemod);
                %%-----------------------2.2, IVA marginal entropy-------------------------------
                if(length(find(wts_blowup))==0) 
                    scv_crosscorr_sqr_sum = zeros(K_scv,1); 
                    scv_ref_crosscorr_sqr_sum = zeros(K_scv,1);  
                    for k = 1:K_scv
                        scv_tmp = [];
                        for dd = 1:length(data)
                           scv_tmp = [scv_tmp,A{dd}(:,k)];  
                        end
                        Sigma{k} = cov(scv_tmp);
                        SCVdevideSigma{k} = scv_tmp/Sigma{k};
                        if find(isnan(SCVdevideSigma{k}))
                            disp(['covariance contains nan for scv', num2str(k)]); 
                        end
                        aSa = mean(sum((SCVdevideSigma{k}).*scv_tmp,2));
                        IVAMarginalEntropy_func = IVAMarginalEntropy_func -log((det(2*pi.*Sigma{k})).^(-1/2))+1/2*aSa;           
                        for dd = 1:length(data)
                            S_inv_deriv_tmp{dd} = 1/(subj_N-1)*(data_org{dd})'*SCVdevideSigma{k}(:,dd);
                            S_inv_deriv_seprt{dd}(:,k) = S_inv_deriv_tmp{dd};                          
                        end

                        [scv_crosscorr_r,scv_crosscorr_p] = corr(scv_tmp,scv_tmp); 
                        scv_crosscorr_r_u_sqr = triu(scv_crosscorr_r,1).^2; 
                        scv_crosscorr_sqr_sum(k) = sum(scv_crosscorr_r_u_sqr(:));
                        scv_ref_crosscorr_sqr_sum(k) = sum(corr(scv_tmp,reference).^2); 
                    end

                    %%%the SCV is selected by sum of the squared correlation between SCV and the reference. 
                    [scv_ref_crosscorr_sqr_sum_max,scv_ref_constrait_ind] = max(scv_ref_crosscorr_sqr_sum);

                    %%%%%%%%%optimize the correlation between the reference and SCV%%%%%%%%%
                    scv_tmp = [];
                    for dd = 1:length(data)
                        scv_tmp = [scv_tmp,A{dd}(:,scv_ref_constrait_ind)]; 
                    end
                    scv_tmp = [scv_tmp,reference];

                    Sigma{scv_ref_constrait_ind} = cov(scv_tmp);
                    SCVdevideSigma{scv_ref_constrait_ind} = scv_tmp/Sigma{scv_ref_constrait_ind};
                    if find(isnan(SCVdevideSigma{scv_ref_constrait_ind}))
                        disp(['covariance contains nan for scv', num2str(scv_ref_constrait_ind)]); 
                    end
                    aSa = mean(sum((SCVdevideSigma{scv_ref_constrait_ind}).*scv_tmp,2));

                    for dd = 1:length(data)
                        S_inv_deriv_tmp{dd} = 1/(subj_N-1)*(data_org{dd})'*SCVdevideSigma{scv_ref_constrait_ind}(:,dd);
                        S_inv_deriv_seprt{dd}(:,scv_ref_constrait_ind) = lambda_ref*S_inv_deriv_tmp{dd};                       
                    end
                end                   
                                  
                %%%%from gradient of S_inv, compute the gradient w.r.t W matrices for each modality
                %%%------------------------------------------------------------------------------------------------------------------------                    
                for ii = 1:length(data)
                    %%%%---------------Gadient for IVA mariginal entropy                   
                    IVAmarginal_diff_weights{ii} = IVAmarginal_diff_weights{ii}-(weights{ii}(p,p)\((S_inv_deriv_seprt{ii})'*VdevideByK{ii}(:,p)))'*weights{ii}(p,p);

                    IVA_diff_weights{ii}(p,p) = -IVA_weights_factor*(IVAjoint_diff_weights{ii}+IVAmarginal_diff_weights{ii});
                    ICA_weight_update_norm = norm(weights_ifmx{ii});
                    fprintf('Data %d, Norm of gradient from ICA %f\n',ii, ICA_weight_update_norm);
                    fprintf('Data %d, Norm of gradient from IVA %f\n',ii, norm(-lambda_IVAref*IVA_diff_weights{ii}));        
                    if (~stopsign(ii))
                        weights{ii} = weights{ii}+ lrate(ii)*(weights_ifmx{ii}+lambda_IVAref*IVA_diff_weights{ii});
                    end
                end
            end                             

            if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i = 1:length(data)
                    weights{i} = weights{i} + momentum*prevwtchange{i};
                    prevwtchange{i} = weights{i}-prevweights{i};
                    prevweights{i} = weights{i};
                end
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if max(max(abs(weights{i}))) > MAX_WEIGHT
                for i = 1:length(data)
                    wts_blowup(i) = 1;
                    change(i) = nochange;
                end
                change_all = nochange;
            end
            if extended 
               for i = 1:length(data)
                if ~wts_blowup(i)
                    %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
                    %
                    if extblocks > 0 & rem(blockno(i),extblocks) == 0,
                        % recompute signs vector using kurtosis
                        if kurtsize(i) < frames(i)
                            rp = randperm(datalength(i));
                            partact=weights{i}*data_pca{1}(:,rp(1:kurtsize(i))); % partial activation
                        else                                                     % for small data sets,
                            partact=weights{i}*data_pca{1};                      % use whole data
                        end
                        m2=mean(partact'.^2).^2;
                        m4= mean(partact'.^4);
                        kk= (m4./m2)-3.0;                           % kurtosis estimates
                        if extmomentum
                            kk = extmomentum*old_kk{i} + (1.0-extmomentum)*kk; % use momentum
                            old_kk{i} = kk;
                        end
                        signs{i}=diag(sign(kk+signsbias));             % pick component signs
                        if signs{i} == oldsigns{i},
                            signcount = signcount+1;
                        else
                            signcount = 0;
                        end
                        oldsigns{i} = signs{i};
                        signcounts = [signcounts signcount];
                        if signcount >= SIGNCOUNT_THRESHOLD,
                            extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                            signcount = 0;                             % less frequent if sign
                        end                                         % is not changing
                    end % extblocks > 0 & . . .                        
                 end  
               end
            end % if extended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for i = 1:length(data)
               blockno(i) = blockno(i) + 1;
            end

            if (length(find(wts_blowup(:)==1))>0)
               break
            end
         end%%%end for all blocks           
         ICAcost_threemod(unique(step)+1) = ICAcost_threemod_sum;
         IVA_cost(unique(step)+1) = -lambda_IVAref*IVA_weights_factor*(IVAMarginalEntropy_func+IVAjoint_func);
    end
    %---------------------------
    % if weight is not  blowup, update
    if length(find(wts_blowup==0)) == length(data)
        delta_concat = [];
        for i = 1:length(data)
           oldwtchange{i} = weights{i}-oldweights{i};
           step(i)=step(i)+1;
           lrates(1,step(i)) = lrate(i);
           angledelta(i)=[0];               
           delta{i}=reshape(oldwtchange{i},1,chans(i)*ncomps(i));
           delta_concat = [delta_concat,delta{i}];
           change(i)=delta{i}*delta{i}';             
        end
           change_all = delta_concat*delta_concat';
     end
    %DATA blow up restart-------------------------------
    if ((length(find(wts_blowup(:)==1))>0) | (~isempty(find(isnan(change)==1))) | (~isempty(find(isinf(change)==1))))% if weights blow up
        fprintf(' ');
        change_all = [0];
        delta_concat = zeros(1,sum(chans.*ncomps));
        olddelta_concat = delta_concat;            
        weights = startweights;
        oldweights = startweights;
        olddelta = delta;
        extblocks = urextblocks;
        prevweights = startweights;

        for i = 1:length(data)
            step(i) = 0;
            stopsign(i)=0;
            if wts_blowup(i) | isnan(change(i))|isinf(change(i))
               lrate(i) = lrate(i)*DEFAULT_RESTART_FAC; % with lower learning rate                 
            end
            change(i) = nochange;
            wts_blowup(i) = 0;    % re-initialize variables
            blockno(i) = [1];
            oldwtchange{i} = zeros(chans(i),ncomps(i));
            delta{i}=zeros(1,chans(i)*ncomps(i));
            prevwtchange{i} = zeros(chans(i),ncomps(i));
            lrates(i,:) = zeros(1,maxsteps);
            bias{i} = zeros(ncomps(i),1);
            if extended
                signs{i} = ones(1,ncomps(i));    % initialize signs to nsub -1, rest +1
                for k=1:nsub
                    signs{i}(k) = -1;
                end
                signs{i} = diag(signs{i}); % make a diagonal matrix
                oldsigns{i} = zeros(size(signs{i}));
            end
            if lrate(i)> MIN_LRATE
                r = rank(data{i});
                if r<ncomps(i)
                    fprintf('For Data %d, Data has rank %d. Cannot compute %d components.\n',...
                        i,r,ncomps(i));
                    return
                else
                    fprintf(...
                        'For Data %d,Lowering learning rate to %g and starting again.\n',i,lrate(i));
                end
            else
                fprintf( ...
                    'runica() : QUITTING - weight matrix may not be invertible! \n');
                return;
            end
        end

    else % if weights in bounds
        for i = 1:length(data)
            %testing the trend of entropy term, avoiding the overfitting of cross-modality correlation optimization               
            lossf1{i}(step(i)) = Compute_entropy_infmx(weights{i},data{i});

            if step(i)>1
                entropychange(i)=lossf1{i}(step(i))-entropy(i);
            else
                entropychange(i)=1;
            end
            
            if step(i)>IVA_weight_adjust_start_step 
                index1(i) = ica_fuse_falsemaxdetect(lossf1{i},trendPara);
            end % end of test------------------------
        end

        if (unique(step)>IVA_weight_adjust_start_step) && (length(find(index1==1))>0)
            disp(['entropy is decreasing for %d modalities ', num2str(length(find(index1==1)))]);
            find(index1==1)
            lambda_IVAref = lambda_IVAref * lambda_IVAref_aneal_fac % anneal weight on IVA with reference if entropy drop  
            if lambda_IVAref < 0.001
                lambda_IVAref = 0.001
            end
            index1 = zeros(1,length(data));
        end

        %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
        if unique(step)> 2 & ~unique(stopsign)
            angledelta_all=acos((delta_concat*olddelta_concat')/sqrt(change_all*oldchange_all));
            fprintf(...
                'Overall,step %d -  wchange %7.6f, angledelta %4.1f deg, IVA_lambda %7.6f, ICA_cost %7.6f, IVA_cost %7.6f\n', ...
                                step(1),change_all,degconst*angledelta_all,lambda_IVAref,ICAcost_threemod(unique(step)),IVA_cost(unique(step)));
        end                            
        for i = 1:length(data)
            if step(i)> 2 & ~stopsign(i)
                angledelta(i)=acos((delta{i}*olddelta{i}')/sqrt(change(i)*oldchange(i)));
            end
            if verbose,
                if step(i) > 2,
                    if ~extended,
                        fprintf(...
                            'Dataset %d, step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
                            i,step(i),lrate(i),change(i),degconst*angledelta(i));
                    else
                        fprintf(...
                            'Dataset %d, step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\n',...
                            i,step(i),lrate(i),change(i),degconst*angledelta(i),(ncomps(i)-sum(diag(signs{i})))/2);
                    end
                elseif ~extended
                    fprintf(...
                        'Dataset %d, step %d - lrate %5f, wchange %7.6f\n',i,step(i),lrate(i),change(i));
                else
                    fprintf(...
                        'Dataset %d, step %d - lrate %5f, wchange %7.6f, %d subgauss\n',...
                        i,step(i),lrate(i),change(i),(ncomps(i)-sum(diag(signs{i})))/2);
                end 
            end; % if verbose
        end  

        %%%%%%%%%%%%%%%%%%%% Anneal learning rate if degree is larger than a certain threshold%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  (~isempty(find(degconst*angledelta > annealdeg))) 
            olddelta  = delta;                % accumulate angledelta until
            oldchange  = change;              % annealdeg is reached
            olddelta_concat = delta_concat;
            oldchange_all = change_all;
            for i =1:length(data)
               if degconst*angledelta(i) > annealdeg
                  lrate(i) = lrate(i)*annealstep; % anneal learning rate
               end
            end
        elseif unique(step) == 1               % on first step only
            olddelta   = delta;                % initialize
            oldchange  = change;
            olddelta_concat = delta_concat;
            oldchange_all = change_all;
        end

        %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (unique(step) >2) & (length(find(change < nochange)) == length(data)) %%% apply stopping rule
            for i = 1:length(data) 
                stopsign(i)=1;   % stop when weights stabilize
            end
        elseif unique(step)>= maxsteps;
            for i = 1:length(data) 
                stopsign(i)=1;   % max step
            end
        elseif (~isempty(find(change > DEFAULT_BLOWUP)))%%% if weights blow up,
            for i =1:length(data)
               if (change(i) > DEFAULT_BLOWUP)
                  lrate(i) = lrate(i)*DEFAULT_BLOWUP_FAC; % with lower learning rate
               end
            end
        end;
        %%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        oldweights = weights;
    end; % end if weights in bounds    

    if (length(find(stopsign==1)) == length(data))
        laststep=step;
        for s = 1:length(step)
            step(s) = maxsteps;
        end                % stop when weights stabilize
        fprintf('aNy-way ICA with reference converged!\n');
     end
%end    
end % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [entropy_y] = Compute_entropy_infmx(w_,x,bias_v)  %%implemented by Jean 
    if exist('bias_v')
        u = w_*x + bias_v*ones(1,size(x,2));
    else
        u = w_*x;
    end
    y=1./(1+exp(-u));
    temp=log(abs(w_*y.*(1-y))+eps);
    entropy_y=mean(temp(:));
end

if (~isempty(find(change > nochange)))
    if max(max(step(:))) == (maxsteps)
        warning('!!!Reached max steps. Please reduce iva_lambda in setup options and restart aNy-way ICA with reference.');
    end
end

if exist('laststep', 'var')
    lrates = lrates(:,1:max(laststep));           % truncate lrate history vector
end

%%%%%%%%%%%%%% Orient components towards positive activation %%%%%%%%%%%
if strcmp(posactflag,'on')
    for i =1:length(data)
        [activations{i},winvout{i},weights{i}] = ica_fuse_posact(data{i},weights{i});
    end
    % changes signs of activations and weights to make activations
    % net rms-positive
else 
    %%%reconstruct final loadings and components
    for i =1:length(data)
        A_f{i} = dewhiteM{i}*pinv(weights{i});
        activations{i} = pinv(A_f{i})*data_org{i};
    end
end
%
%%%%%%%%%%%%%% If pcaflag, compose PCA and ICA matrices %%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Composing the eigenvector, weights, and sphere matrices\n');
    fprintf('  into a single rectangular weights matrix; sphere=eye(%d)\n'...
        ,chans);
    for i =1:length(data)
        weights{i}= weights{i}*sphere{i}*eigenvectors{i}(:,1:ncomps(i))';
        sphere{i} = eye(urchans(i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end
