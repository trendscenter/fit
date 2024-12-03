function [weightsF,IVAref_stepF,IVAref_cost] = ica_fuse_IVAref_Nway(data,data_org,reference,VdevideByK,lambda_adjust_fac,dewhiteM)
%%%data: dimension of PC-by-variable(prewhitened) 
%%%data_org: dimension of subject-by-variable
%%%reference: dimension of subject-by-one
%%%VdevideByK: obtained from SVD decomposition when doing dimension reduction

IC_num = [];
for m = 1:length(data)
    [chans(m) frames(m)] = size(data{m}); % determine the data size
    IC_num = [IC_num,size(data{m},1)];
end
IC_min = min(IC_num);
subj_N = size(data_org{1},1);

%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      % anneal by multiplying lrate by this original 0,9 to 0.95 changed by JL

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       %%DEFAULT_BLOWUP_FAC   = 0.8; % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.95;       % %%DEFAULT_RESTART_FAC  = 0.9;  if weights blowup, restart with lrate
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_LRATE        = 0.015./log(chans); 
DEFAULT_EXTENDED     = 0;         % default off

%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
extended   = DEFAULT_EXTENDED;
verbose    = DEFAULT_VERBOSE;
lrate      = DEFAULT_LRATE;
annealdeg  = DEFAULT_ANNEALDEG;
annealstep = DEFAULT_ANNEALSTEP;     % defaults defined above
nochange   = DEFAULT_STOP;
ncomps     = chans;
wts_blowup = zeros(1,length(data));  % flag =1 when weights too large
maxsteps = 512;  %%%%training after this many steps

weights=[];delta = []; prevwtchange = []; oldwtchange = [];
for i = 1:length(data)
   weights{i} = eye(ncomps(i),chans(i)); % begin with the identity matrix
   delta{i}=zeros(1,chans(i)*ncomps(i));
   prevwtchange{i} = zeros(chans(i),ncomps(i));
   oldwtchange{i} = zeros(chans(i),ncomps(i));
end
weightsF = weights;
IVAref_stepF = 0;
degconst = 180./pi;
startweights = weights;
oldweights = startweights;
delta_concat = zeros(1,sum(chans.*ncomps));

step_IVAref=0;
stopsign_IVAref=zeros(1,length(data));      
angledelta=zeros(1,length(data));
change=zeros(1,length(data));
                
p = 1:min(ncomps);
K_scv = min(ncomps); 
lambda_IVAref_preprocess = lambda_adjust_fac*300*max(lrate);
lambda_ref = 0.001;
%%%%%%%-------------------------------------------------------------------------------------------
%%%%%%%Apply IVA with reference to pre-align the SCVs to provide a good starting point for aNy-way ICA with reference
while step_IVAref < maxsteps
    if ~isempty(find(stopsign_IVAref==0))
        for i = 1:length(data)
           SInv{i} = VdevideByK{i}/oldweights{i}; 
        end
        IVAref_jointFunc = 0;
        for i = 1:length(data)           
            B{i} = SInv{i}(:,p);
            [V2{i},S2{i}] = eig(B{i}'*B{i});
            S2{i} = sqrt(diag(S2{i}));
            U2{i} = B{i}/((V2{i})'.*repmat(S2{i},1,size(V2{i},2)));
            IVAref_jointFunc = IVAref_jointFunc-sum(log(abs(S2{i}(:))));
            Q{i} = zeros(K_scv,K_scv);
            Q{i}(p,:) = V2{i};

            %-------------------------------------------
            F{i} = U2{i}'*VdevideByK{i}(:,p); 
            FdevideW{i} = F{i}/oldweights{i}(p,p);
            T{i} = FdevideW{i}*Q{i};
            QdevideT{i} = Q{i}/T{i}; 
            %%%%-----------natural gradient(regular gradient*W^T*W)
            IVArefJoint_diff_weights{i} = (oldweights{i}(p,p)\((QdevideT{i})*F{i}))'*oldweights{i}(p,p); 
        end

        IVArefMarginalEntropy_func = 0;
        for i = 1:length(data)
           IVAref_marginal_diff_weights{i} = zeros(K_scv,K_scv);
           A{i} = data_org{i}*SInv{i};
           S_inv_deriv_seprt{i} = zeros(frames(i),K_scv);
        end
        
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
            IVArefMarginalEntropy_func = IVArefMarginalEntropy_func -log((det(2*pi.*Sigma{k})).^(-1/2))+1/2*aSa;           
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

        %%%%%%%%%optimize the correlation between the reference and SCV
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
       
        IVAref_cost(step_IVAref+1) = -lambda_IVAref_preprocess*(IVArefMarginalEntropy_func+IVAref_jointFunc);

        for i = 1:length(data)
            if (~stopsign_IVAref(i))
                %%%%natural gradient(regular gradient*W^T*W)
                IVAref_marginal_diff_weights{i} = IVAref_marginal_diff_weights{i}-(oldweights{i}(p,p)\((S_inv_deriv_seprt{i})'*VdevideByK{i}(:,p)))'*oldweights{i}(p,p);
                IVAref_diff_weights{i} = -(IVArefJoint_diff_weights{i}+IVAref_marginal_diff_weights{i});
                weights{i}(p,p) = weights{i}(p,p)+ lambda_IVAref_preprocess*IVAref_diff_weights{i};
            end
            if max(max(abs(weights{i}))) > MAX_WEIGHT 
                wts_blowup(i) = 1;
                change(i) = nochange;
                change_all = 0;
            end
        end

        %---------------------------
        % if weight is not  blowup, update
        if length(find(wts_blowup==0)) == length(data)
            delta_concat = [];
            step_IVAref = step_IVAref+1;
            for i = 1:length(data)
                oldwtchange{i} = weights{i}-oldweights{i};               
                angledelta(i)=[0];               
                delta{i}=reshape(oldwtchange{i},1,chans(i)*ncomps(i));
                delta_concat = [delta_concat,delta{i}];
                change(i)=delta{i}*delta{i}';             
            end
            change_all = delta_concat*delta_concat';

            if step_IVAref>50 
                if change_all<=change_all_prev 
                    weightsF = weights;  
                    IVAref_stepF = step_IVAref;
                    change_all_prev = change_all;
                else
                    if step_IVAref>IVAref_stepF+10 
                        step_IVAref = maxsteps;
                        for ii = 1:length(data)
                            stopsign_IVAref(ii) = 1;
                        end
                    end
                end
            else
                weightsF = weights;
                IVAref_stepF = step_IVAref;
                change_all_prev = change_all;
            end
        end

        if ((length(find(wts_blowup(:)==1))>0) | (~isempty(find(isnan(change)==1))) | (~isempty(find(isinf(change)==1))))% if weights blow up,
            fprintf(' ');
            change_all = [0];
            change_all_prev = [0];
            delta_concat = zeros(1,sum(chans.*ncomps));
            olddelta_concat = delta_concat;            
            weights = startweights;        % and original weight matrix
            oldweights = startweights;
            olddelta = delta;
            lambda_IVAref_preprocess = lambda_IVAref_preprocess*DEFAULT_RESTART_FAC; % with lower learning rate  
            step_IVAref = 0;
            fprintf('Lowering lambda to %g and starting again.\n',lambda_IVAref_preprocess);
            for i = 1:length(data)                
                stopsign_IVAref(i)=0;
                change(i) = nochange;
                wts_blowup(i) = 0;    % re-initialize variables
                oldwtchange{i} = zeros(chans(i),ncomps(i));
                delta{i}=zeros(1,chans(i)*ncomps(i));
                prevwtchange{i} = zeros(chans(i),ncomps(i));
                if lrate(i)> MIN_LRATE
                    r = rank(data{i});
                    if r<ncomps(i)
                        fprintf('For Data %d, Data has rank %d. Cannot compute %d components.\n',...
                            i,r,ncomps(i));
                        return
                    end
                else
                    fprintf( ...
                        'runica() : QUITTING - weight matrix may not be invertible! \n');
                    return;
                end
            end
            
        else      
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step_IVAref > 2 & ~unique(stopsign_IVAref)
                angledelta_all=acos((delta_concat*olddelta_concat')/sqrt(change_all*oldchange_all));
                fprintf(...
                    'Overall,step_IVAref %d -  wchange %7.6f, angledelta %4.1f deg\n', ...
                                    step_IVAref,change_all,degconst*angledelta_all);
            end                            
            for i = 1:length(data)
                if step_IVAref > 2 & ~stopsign_IVAref(i)
                    angledelta(i)=acos((delta{i}*olddelta{i}')/sqrt(change(i)*oldchange(i)));
                end
                if verbose,
                    if step_IVAref > 2,
                        if ~extended,
                            fprintf(...
                                'Dataset %d, step_IVAref %d - IVA_lambda %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
                                i,step_IVAref,lambda_IVAref_preprocess,change(i),degconst*angledelta(i));
                        else
                            fprintf(...
                                'Dataset %d, step_IVAref %d - IVA_lambda %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\n',...
                                i,step_IVAref,lambda_IVAref_preprocess,change(i),degconst*angledelta(i),(ncomps(i)-sum(diag(signs{i})))/2);
                        end
                    elseif ~extended
                        fprintf(...
                            'Dataset %d, step_IVAref %d - IVA_lambda %5f, wchange %7.6f\n',i,step_IVAref,lrate(i),change(i));
                    else
                        fprintf(...
                            'Dataset %d, step_IVAref %d - IVA_lambda %5f, wchange %7.6f, %d subgauss\n',...
                            i,step_IVAref,lambda_IVAref_preprocess,change(i),(ncomps(i)-sum(diag(signs{i})))/2);
                    end 
                end; 
            end
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate if degree is larger than the threshold%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if  (~isempty(find(degconst*angledelta > annealdeg)))
                olddelta  = delta;                % accumulate angledelta until
                oldchange  = change;               %  annealdeg is reached
                olddelta_concat = delta_concat;
                oldchange_all = change_all;
                lambda_IVAref_preprocess = lambda_IVAref_preprocess*annealstep; % anneal learning rate                
            elseif step_IVAref == 1                     
                olddelta   = delta;                % initialize
                oldchange  = change;
                olddelta_concat = delta_concat;
                oldchange_all = change_all;
            end
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (step_IVAref >2) & (length(find(change < nochange)) == length(data))  % apply stopping rule
                for i = 1:length(data) 
                    stopsign_IVAref(i)=1;   % stop when weights stabilize
                end
            elseif step_IVAref >= maxsteps;
                for i = 1:length(data) 
                    stopsign_IVAref(i)=1;   % stop when max step_IVAref
                end               
            elseif (~isempty(find(change > DEFAULT_BLOWUP)))    % if weights blow up,
                lambda_IVAref_preprocess = lambda_IVAref_preprocess*DEFAULT_BLOWUP_FAC; % with lower learning rate
            end  
            %%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            oldweights = weights;  
        end    
    end
    if (length(find(stopsign_IVAref==1)) == length(data))
        step_IVAref = maxsteps;
        fprintf('IVA with reference converged!\n');
    end
end



              
         