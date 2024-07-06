% 082523, Cyrus Eierud
% Below is a new algorithm developed by Ibrahim Khalilullah
% the new algo may use parts of Infomax, but also uses machine learning
%
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
% 'lrate'     = [rate] initial ICA learning rate (<< 1) (default -> heuristic)
% 'anneal'    = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended)
%                         controls speed of convergence
% 'stop'      = [f] stop training when weight-change < this (default -> 1e-6)
% 'maxsteps'  = [N] max number of ICA training steps    (default -> 512)
% 'bias'      = ['on'/'off'] perform bias adjustment    (default -> 'on')
% 'momentum'  = [0<f<1] training momentum               (default -> 0)
% 'extended'  = [N] perform tanh() "extended-ICA" with sign estimation 
%               every N training blocks. If N < 0, fix number of sub-Gaussian
%               components to -N [faster than N>0]      (default|0 -> off)
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

function [weights,sphereGmCm]=ica_fuse_OptInfomax(dataSeparate,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)
    if nargin < 1
       help runica  
       return
    end
    
    [chans frames] = size(dataSeparate); % determine the data size
    urchans = chans;  % remember original data channels 
    datalength = frames;
    %
    %%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
    %
    DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
    DEFAULT_ANNEALSTEP   = 0.98;      %     anneal by multiplying lrate by this
    DEFAULT_MAXSTEPS     = 100;       % ]top training after this many steps 
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
    
    DEFAULT_SPHEREFLAG   = 'on';      % use the sphere matrix as the default
    %   starting weight matrix
    DEFAULT_PCAFLAG      = 'off';     % don't use PCA reduction
    DEFAULT_POSACTFLAG   = 'on';      % use posact()
    DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
    DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule
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
    annealstep = 0;                      % defaults declared below
    nochange   = DEFAULT_STOP;
    maxsteps   = DEFAULT_MAXSTEPS;
    
    ncomps     = chans;
    
    extblocks  = DEFAULT_EXTBLOCKS;
    kurtsize   = MAX_KURTSIZE;
    signsbias  = 0.02;                   % bias towards super-Gaussian components
    extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
    nsub       = DEFAULT_NSUB;
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
       
       if strcmp(Keyword,'ncomps')
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
       elseif strcmpi(Keyword,'sPcaFile')
          if ~isstr(Value)
             fprintf('runica(): sPcaFile must be a string')
             return
          end
          sPcaFile = Value;
       elseif strcmpi(Keyword,'iMod1LenOnes')
          if isstr(Value)
             fprintf('runica(): iMod1LenOnes must be a Value')
             return
          end
          iMod1LenOnes = Value;            
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
       else
          fprintf('runica(): unknown flag')
          return
       end
    end
    %
    
    dataSeparate = dataSeparate';
    num_mods = 2; % Supports only 2 modalities so  far
    ann_stp = annealstep;
    max_steps=maxsteps;
    clear maxsteps;
    
    whitesigGm=dataSeparate(1:iMod1LenOnes, :);  
    whitesigCm=dataSeparate(iMod1LenOnes+1:end,:); 
    
    dataGm=whitesigGm';
    [chansGm, framesGm] = size(dataGm); % determine the data size
    dataCm=whitesigCm';
    [chansCm, framesCm] = size(dataCm); % determine the data size
    urchansCm = chansCm;
    
    
    %%% initialize weights
    weights = randsmall(ncomps,chansGm);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pass=1;    
    
    %%%%%%% imortant for termination and optimization %%%%%%%%
    degconst = 180./pi;
    annealstep = 0.97;
    annealdeg  = 60;
    %%%%%%% imortant for termination and optimization %%%%%%%%       
    
    DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
    DEFAULT_BLOWUP_FAC   = 0.9;
    posactflag='on';
    wts_passed=0;
    
    figGauge = figure,
    
    while pass<=max_steps
        if pass==1
            maxsteps=1;
        elseif pass>1 && pass<10
            maxsteps=2;
        elseif pass>=10 && pass<20
            maxsteps=5;
        elseif pass>=20 && pass<30
            maxsteps=10;
        elseif pass>=30 && pass<40
            maxsteps=15;
        elseif pass>=40 && pass<50
            maxsteps=20;
        elseif pass>=50
            maxsteps=maxsteps+2;
        end
        
        [W1, lrate, ~] = ica_fuse_pmljICA(dataGm, ncomps, weights, chansGm, framesGm, lrate, maxsteps);
        weights=W1; 
        
        [W2, lrate, sphere, data, signs, bias] = ica_fuse_pmljICA(dataCm, ncomps, weights, chansCm, framesCm, lrate, maxsteps);

        weights=(W1+W2)./num_mods;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pass==1
            oldweights=weights;   
            oldwtchange=W2-W1;
            olddelta=reshape(oldwtchange,1,chansCm*ncomps);
            oldchange = olddelta*olddelta';   
            change=999;
    
        else
    
            oldwtchange=weights-oldweights;
            delta=reshape(oldwtchange,1,chansCm*ncomps);
            change=delta*delta'; 
            angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
            %
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if degconst*angledelta > annealdeg, 
               if lrate<0.00001
                   annealstep = ann_stp;%%               
               end
               
               lrate = lrate*annealstep;          % anneal learning rate
               olddelta   = delta;                % accumulate angledelta until
               oldchange  = change;               %  annealdeg is reached                      
            end
    
            fprintf(...
                      'step %d - lrate %7.9f, wchange %7.7f, angledelta %4.1f deg\n',...
                      pass,lrate,change,degconst*angledelta);
    
            oldweights=weights;
            
            if verbose==1
                if pass==2
                    plt_change(pass)=change;
                    plot(2:pass, plt_change(2:pass), '-ob');
                    set(gca,'XTick',(1:50:701))
                    grid on
                    title('Jointly weights optimization')
                    xlabel('Epoch')
                    ylabel('Weights')
                    hold on
                elseif pass<max_steps
                    plt_change(pass)=change;
                    plot(2:pass, plt_change(2:pass), '-ob');   
                    set(gca,'XTick',(1:50:701))
                    pause(0.01)
                    grid on
                    title('Joint Optimization')
                    xlabel('Epoch')
                    ylabel('Weight change')
                end
            end
    
        end
        
        pass=pass+1;    
        
        if pass >2 & change < nochange,      % apply stopping rule
           laststep=pass;            
           pass=max_steps+1;                  % stop when weights stabilize
        elseif change > DEFAULT_BLOWUP,      % if weights blow up,
           lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying 
        end; 
    end
    close(figGauge);
    if verbose==1
        hold off;
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu1=mean(dataSeparate);  
    
    type_pca='standard';
    load(sPcaFile, 'whitesigCombined', 'dewhiteM');
    whitesigGmCm = whitesigCombined';
    dewhiteMGmCm = dewhiteM';
    
    dataGmCm=whitesigGmCm';
    
    %% comment it if you sphereing before or use PCA sphering
    sphereGmCm = 2.0*inv(sqrtm(cov(dataGmCm'))); 
    
    %%%%%%%%%%%%%% Orient components towards positive activation %%%%%%%%%%%
    data = sphereGmCm*dataGmCm; 
    if strcmp(posactflag,'on')
       [activations, winvout, weights] = ica_fuse_posact(data, weights);
       % changes signs of activations and weights to make activations
       % net rms-positive
    else
       activations = weights*data;
    end
    
    fprintf(...
          'Sorting components in descending order of mean projected variance ...\n');  
    
    if wts_passed == 0
       %
       %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %
       meanvar  = zeros(ncomps,1);      % size of the projections
       if ncomps == urchansCm % if weights are square . . . %% equal to number of component
          winv = inv(weights*sphereGmCm);
       else
          fprintf('Using pseudo-inverse of weight matrix to rank order component projections.\n');
          winv = pinv(weights*sphereGmCm);
       end
       for s=1:ncomps
          
          fprintf('%d ',s);         % construct single-component data matrix      
          % project to scalp, then add row means load ica_fusion_oasis_53ICN_avgwtsh_30ica2_std.mat;
          compproj = winv(:,s)*activations(s,:);
          meanvar(s) = mean(sum(compproj.*compproj)/(size(compproj,1)-1));
          % compute mean variance 
       end                                         % at all scalp channels
       
       fprintf('\n');
       %
       %%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
       %
       [sortvar, windex] = sort(meanvar);
       windex = windex(ncomps:-1:1); % order large to small 
       meanvar = meanvar(windex);
       % 
       %%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
       %
       
       fprintf('Permuting the activation wave forms ...\n');   
       activations = activations(windex,:);
       
       weights = weights(windex,:);% reorder the weight matrix
       
    else
       fprintf('Components not ordered by variance.\n');
    end
    %
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     W1=weights*sphereGmCm;
%     sources=W1*dataGmCm;
%     A = pinv(W1);
%     
%     A = dewhiteMGmCm*A;
%     W = pinv(A);
    
end


