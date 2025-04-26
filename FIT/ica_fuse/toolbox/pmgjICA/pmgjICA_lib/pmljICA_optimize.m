

%%%function [weights, sphere, bias, signs, lrates, activations]=my_jICA(data, ncomps, weights, chans, frames)
function [weights, lrate, sphere, data, signs, bias]=pmljICA_optimize(data, maxsteps, lrate, weights, p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)

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
%%%%DEFAULT_MAXSTEPS     = 512;       % ]top training after this many steps 
DEFAULT_MOMENTUM     = 0.0;       % default momentum weight

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.9;       % if weights blowup, restart with lrate
% lower by this factor
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
%%DEFAULT_LRATE        = 0.015/log(chans); 
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
%%%%lrate      = DEFAULT_LRATE;
annealdeg  = DEFAULT_ANNEALDEG;
annealstep = 0;                      % defaults declared below
nochange   = DEFAULT_STOP;
momentum   = DEFAULT_MOMENTUM;
%%%%maxsteps   = DEFAULT_MAXSTEPS;

%%weights    = 0;                      % defaults defined below
ncomps     = chans;
biasflag   = DEFAULT_BIASFLAG;

extended   = DEFAULT_EXTENDED;
extblocks  = DEFAULT_EXTBLOCKS;
kurtsize   = MAX_KURTSIZE;
signsbias  = 0.02;                   % bias towards super-Gaussian components
extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
nsub       = DEFAULT_NSUB;
wts_blowup = 0;                      % flag =1 when weights too large
wts_passed = 0;                      % flag weights passed as argument
%
%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
% if (nargin> 1 & rem(nargin,2) == 0)
%    fprintf('runica(): Even number of input arguments???')
%    return
% end
for i = 6:2:nargin % for each Keyword
   Keyword = eval(['p',int2str((i-6)/2 +1)]);
   Value = eval(['v',int2str((i-6)/2 +1)]);
   if ~isstr(Keyword)
      fprintf('runica(): keywords must be strings')
      return
   end
   
      
   if strcmp(Keyword,'ncomps')
      if isstr(Value)
         fprintf('runica(): ncomps value must be an integer')
         return
      end
      if ncomps < urchans & ncomps ~= Value
         fprintf('runica(): Use either PCA or ICA dimension reduction');
         return
      end
      
      if strcmp(pcaflag, 'off') & ncomps
          ncomps=chans;
      else
          ncomps = Value;
      end

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
      pcaflag = 'off';

      if strcmp(pcaflag, 'off') & ~ncomps
          ncomps = Value;
      end
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
       %%%%disp('Inside verbose    ..... ')
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
      verbose=0;  %%% always zero
   %%else
      %%fprintf('runica(): unknown flag')
      %%return
   end
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
%%if verbose,
%%   fprintf('Removing mean of each channel ...\n');
%%%end
%%% disp('Not removing mean of each channel!!!');
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
   [eigenvectors,eigenvalues,data] = mygift_pcsquash(data,ncomps);
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
   if verbose,
      fprintf('Sphering the data ...\n');
   end
   data = sphere*data;      % actually decorrelate the electrode signals
   
elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  if verbose,
     fprintf('Using starting weights named on commandline ...\n');
     fprintf('Returning the identity matrix in variable "sphere" ...\n');
  end
  sphere = eye(chans);                 % return the identity matrix
   
elseif strcmp(sphering,'none')
 
  if verbose,
     fprintf('Using starting weights named on commandline ...\n');
     fprintf('Returning the identity matrix in variable "sphere" ...\n');
  end
   
  sphere = eye(chans,chans);
  if verbose,
     fprintf('Returned variable "sphere" will be the identity matrix.\n');
  end
end

%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%

lastt=fix((datalength/block-1)*block+1);
BI=block*eye(ncomps,ncomps);
delta=zeros(1,chans*ncomps);
changes = [];
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




%disp('block.....')
%%%%%divide
%block
%weights


while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   permute=randperm(datalength); % shuffle data order at each step
     
   %%disp('weights...')
   %%%size(weights)

   for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
      
      if biasflag                                                   
         u=weights*data(:,permute(t:t+block-1)) + bias*onesrow;      
      else                                                             
         u=weights*data(:,permute(t:t+block-1));                      
      end   
      
      %trainblk=size(u)  %%% 5 595
      
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

      blockno = blockno + 1;
      if wts_blowup
         break
      end
   end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if ~wts_blowup
      oldwtchange = weights-oldweights;
      step=step+1; 
      %
      %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
      %
      lrates(1,step) = lrate;
      angledelta=0.;
      delta=reshape(oldwtchange,1,chans*ncomps);
      change=delta*delta'; 
   end

   %
   %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
   %
   if wts_blowup | isnan(change)|isinf(change),  % if weights blow up,
      fprintf('... Weights blowup!!!');
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
         %%lrate
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
      end
      
      oldweights = weights;
      %
      %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      if degconst*angledelta > annealdeg,  
         if lrate<0.000001
               annealstep = 0.96;%%
         end
         lrate = lrate*annealstep;          % anneal learning rate
         olddelta   = delta;                % accumulate angledelta until
         oldchange  = change;               %  annealdeg is reached
      elseif step == 1                     % on first step only
         olddelta   = delta;                % initialize 
         oldchange  = change;               
      end
      %
      %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      if step >2 & change < nochange,      % apply stopping rule                    
         step=maxsteps;                  % stop when weights stabilize
      elseif change > DEFAULT_BLOWUP,      % if weights blow up,
         lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying 
      end;

   end
   
end; % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


