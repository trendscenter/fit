

%%%function [weights, sphere, bias, signs, lrates, activations]=my_jICA(data, ncomps, weights, chans, frames)
function [weights, lrate, sphere, data, signs, bias]=my_jICA(data, ncomps, weights, chans, frames, lrate, maxsteps)

degconst = 180./pi;
datalength = frames;
DEFAULT_BLOWUP_in       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC_in   = 0.9;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC_in  = 0.9;       % if weights blowup, restart with lrate
annealdeg_in = 60;

%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
annealstep_in = 0.97;
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up

% lower by this factor
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit

block      = floor(sqrt(frames/3));          % heuristic default - may need adjustment!
                      % defaults declared below
nochange   = 0.000001;  % stop training if weight changes below this
momentum   = 0.0;       % default momentum weight
biasflag   = 1;
extended   = 0;
extblocks  = 1;
nsub       = 1;
wts_blowup = 0;                      % flag =1 when weights too large


sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
data = sphere*data;      % actually decorrelate the electrode signals


lastt=fix((datalength/block-1)*block+1);
BI=block*eye(ncomps,ncomps);
startweights = weights;
prevweights = startweights;
oldweights = startweights;
prevwtchange = zeros(chans,ncomps);

lrates = zeros(1,maxsteps);
onesrow = ones(1,block);
bias = zeros(ncomps,1);
signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1

for k=1:nsub
   signs(k) = -1;
end

signs = diag(signs); % make a diagonal matrix
% % % oldsigns = zeros(size(signs));
% % % signcount = 0;              % counter for same-signs
% % % signcounts = [];
urextblocks = extblocks;    % original value, for resets
% % % % old_kk = zeros(1,ncomps);   % for kurtosis momemtum

step=0;
blockno = 1;  % running block counter for kurtosis interrupts



%%%%%divide

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
      lrate = lrate*DEFAULT_RESTART_FAC_in; % with lower learning rate
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
      if degconst*angledelta > annealdeg_in,  
         if lrate<0.00001
               annealstep_in = 0.98;%%
         end
         lrate = lrate*annealstep_in;          % anneal learning rate
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
      elseif change > DEFAULT_BLOWUP_in,      % if weights blow up,
         lrate=lrate*DEFAULT_BLOWUP_FAC_in;    % keep trying 
      end;

   end
   
end; % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


