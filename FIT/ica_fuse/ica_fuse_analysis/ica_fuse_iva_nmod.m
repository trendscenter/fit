function [weightsF,iva_stepF,IVA_cost,weight_change_overall] = ica_fuse_iva_nmod(data,data_org,VdevideByK,lambda_adjust_fac)
%%%data: dimension of PC-by-variable(whitened data) 
%%%data_org: dimension of subject-by-variable
%%%VdevideByK: obtained from SVD decomposition when doing dimension
%%%reduction
%--------------------------------------------------------------------------
for i = 1:length(data)
    [chans(i) frames(i)] = size(data{i}); % determine the data size
end
subj_N = size(data_org{1},1);

%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      % anneal by multiplying lrate by this original 0,9 to 0.95 changed by JL
DEFAULT_MAXSTEPS     = 512;       %training after this many steps 512

DEFAULT_BLOWUP       = 1000000000.0;  % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.95;      % if weights blowup, restart with lrate
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
% maxsteps   = DEFAULT_MAXSTEPS;
ncomps     = chans;
maxsteps = 512;
wts_blowup = zeros(1,length(data));   % flag =1 when weights too large

%%%%%%%%%%%%%%%%%%%%%%%% Initialize IVA training %%%%%%%%%%%%%%%%%%%%%%%%%
weights=[];
delta = [];
prevwtchange = [];
oldwtchange = [];
for i = 1:length(data)
    weights{i} = eye(ncomps(i),chans(i)); % begin with the identity matrix
    delta{i}=zeros(1,chans(i)*ncomps(i));
    prevwtchange{i} = zeros(chans(i),ncomps(i));
    oldwtchange{i} = zeros(chans(i),ncomps(i));
end
weightsF = weights;
startweights = weights;
oldweights = startweights;

stopsign_IVA = zeros(1,length(data));   % flag =1 when it converges
angledelta= zeros(1,length(data));
change = zeros(1,length(data));

degconst = 180./pi;
delta_concat = zeros(1,sum(chans.*ncomps));

iva_stepF = 0;
step_IVA=0;                
%%%%%%%%set minimum component number and component indices for IVA to
%%%%%%%%generate SCV
p = 1:min(ncomps);
K_scv = min(ncomps); 
lambda_IVA_preprocess = lambda_adjust_fac*300*max(lrate);

%-------------------------------------------------------------------------------------------------
%%%%%%%Apply IVA only to pre-align the SCVs to provide a good starting point
%%%%%%%for aNy-way ICA
while step_IVA < maxsteps
    if ~isempty(find(stopsign_IVA==0))   %%%%(~stopsign_IVA(1))|(~stopsign_IVA(2))|(~stopsign_IVA(3))
        %%%compute inverse of S based on weights from last step_IVA for IVA gradient calculation
        for m = 1:length(data)
            SInv{m} = VdevideByK{m}/oldweights{m};%%%Update inverse of S with updated weight matrix
        end

        IVAjoint_func = 0;
        %%%derivation for joint entropy part
        for i = 1:length(data)           
            B{i} = SInv{i}(:,p);%%%pick components for correlation optimization.
            [V2{i},S2{i}] = eig(B{i}'*B{i});
            S2{i} = sqrt(diag(S2{i}));
            U2{i} = B{i}/((V2{i})'.*repmat(S2{i},1,size(V2{i},2)));
            IVAjoint_func = IVAjoint_func-sum(log(abs(S2{i}(:))));
            Q{i} = zeros(K_scv,K_scv);
            Q{i}(p,:) = V2{i};%%extend from common component number dimension to original component dimendion for each modilaty

            %-------------------------------------------
            F{i} = U2{i}'*VdevideByK{i}(:,p); 
            FdevideW{i} = F{i}/oldweights{i}(p,p);
            T{i} = FdevideW{i}*Q{i};
            QdevideT{i} = Q{i}/T{i}; 

            %%%%-----------natural gradient(regular gradient*W^T*W)
            IVAjoint_diff_weights{i} = (oldweights{i}(p,p)\((QdevideT{i})*F{i}))'*oldweights{i}(p,p); 
        end
        %%%%derivation for marginal entropy part 
        
        for i = 1:length(data) 
            IVAmarginal_diff_weights{i} = zeros(K_scv,K_scv);
            A{i} = data_org{i}*SInv{i};%%A{i} = dewhiteM{i}*pinv(oldweights{i}); %%% 
            S_inv_deriv_seprt{i} = zeros(frames(i),K_scv);
        end
        
        IVAMarginalEntropy_func = 0;
        for k = 1:K_scv
            scv_tmp = [];
            for i = 1:length(data)
                scv_tmp = [scv_tmp,A{i}(:,k)];
            end
            Sigma{k} = cov(scv_tmp);
            SCVdevideSigma{k} = scv_tmp/Sigma{k};
            if find(isnan(SCVdevideSigma{k}))
                disp(['covariance contains nan for scv', num2str(k)]); 
            end
            aSa = mean(sum((SCVdevideSigma{k}).*scv_tmp,2));

            IVAMarginalEntropy_func = IVAMarginalEntropy_func -log((det(2*pi.*Sigma{k})).^(-1/2))+1/2*aSa; 
            
            for i = 1:length(data)
                S_inv_deriv_tmp{i} = 1/(subj_N-1)*(data_org{i})'*SCVdevideSigma{k}(:,i);
                S_inv_deriv_seprt{i}(:,k) = S_inv_deriv_tmp{i};
            end
        end

        IVA_cost(step_IVA+1) = -lambda_IVA_preprocess*(IVAMarginalEntropy_func+IVAjoint_func);

        for i = 1:length(data)
            if (~stopsign_IVA(i))
                %%%---------------Gadient for IVA mariginal entropy-----------
                %%%%-----------natural gradient(regular gradient*W^T*W)-----------
                IVAmarginal_diff_weights{i} = IVAmarginal_diff_weights{i}-(oldweights{i}(p,p)\((S_inv_deriv_seprt{i})'*VdevideByK{i}(:,p)))'*oldweights{i}(p,p);
                %%%%---------------Gadient for IVA mariginal entropy

                %%%%%%-----------Gradient update for IVA part-----------
                IVA_diff_weights{i} = -(IVAjoint_diff_weights{i}+IVAmarginal_diff_weights{i});
                %%%%%%-----------Gradient update for IVA part-----------
                weights{i}(p,p) = weights{i}(p,p)+ lambda_IVA_preprocess*IVA_diff_weights{i};
            end
            if max(max(abs(weights{i}))) > MAX_WEIGHT 
                wts_blowup(i) = 1;
                change(i) = nochange;
                change_all = nochange;
            end
        end

        %---------------------------
        % if weight does not  blowup, update
        if length(find(wts_blowup==0)) == length(data) %%%(~wts_blowup(1))&(~wts_blowup(2))&(~wts_blowup(3))
            delta_concat = [];
            step_IVA = step_IVA+1;
            for i = 1:length(data)
                oldwtchange{i} = weights{i}-oldweights{i};               
                angledelta(i)= 0;               
                delta{i}=reshape(oldwtchange{i},1,chans(i)*ncomps(i));
                delta_concat = [delta_concat,delta{i}];
                change(i)=delta{i}*delta{i}';             
            end
            change_all = delta_concat*delta_concat';
            weight_change_overall(step_IVA+1) = change_all;
            
            %%%%IVA find the weight with the minimum weight update
            if step_IVA>50 %%%%let IVA run for 50 steps and allow for vibration
                if change_all<=change_all_prev %%%if the weight update keeps decreasing, assign the weightF with the current weight,otherwise, keep monitoring weight change for 10 more steps, if it continues increasing, stop IVA iterations
                    weightsF = weights;  
                    iva_stepF = step_IVA;
                    change_all_prev = change_all;
                else
                    if step_IVA>iva_stepF+10 %%%%
                        step_IVA = maxsteps;
                        for i = 1:length(data)
                            stopsign_IVA(i) = 1;
                        end
                    end
                end
            else
                weightsF = weights;
                iva_stepF = step_IVA;
                change_all_prev = change_all;
            end
        end
        
        if ((length(find(wts_blowup(:)==1))>0) | (~isempty(find(isnan(change)==1))) | (~isempty(find(isinf(change)==1)))) %%%((length(find(wts_blowup(:)==1))>0) | (isnan(change(1))|isnan(change(2))|isnan(change(3)))|(isinf(change(1))|isinf(change(2))|isinf(change(3)))) % if weights blow up,
            fprintf(' ');
            change_all = [0];
            change_all_prev = [0];   
            delta_concat = zeros(1,sum(chans.*ncomps));                        
            olddelta_concat = delta_concat;            
            weights = startweights;    
            oldweights = startweights;
            olddelta = delta;
            lambda_IVA_preprocess = lambda_IVA_preprocess*DEFAULT_RESTART_FAC; % with lower learning rate  
            step_IVA = 0;
            fprintf('Lowering lambda to %g and starting again.\n',lambda_IVA_preprocess);
            for i = 1:length(data)                
                stopsign_IVA(i)=0;
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
            
        else % if DATA1 weights in bounds           
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            if step_IVA > 2 & ~unique(stopsign_IVA)
                angledelta_all=acos((delta_concat*olddelta_concat')/sqrt(change_all*oldchange_all));
                fprintf(...
                    'Overall,step_IVA %d -  wchange %7.6f, angledelta %4.1f deg\n', ...
                                    step_IVA,change_all,degconst*angledelta_all);
            end                            
            for i = 1:length(data)
                if step_IVA > 2 & ~stopsign_IVA(i)
                    angledelta(i)=acos((delta{i}*olddelta{i}')/sqrt(change(i)*oldchange(i)));
                end
                if verbose,
                    if step_IVA > 2,
                        if ~extended,
                            fprintf(...
                                'Dataset %d, step_IVA %d - IVA_lambda %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
                                i,step_IVA,lambda_IVA_preprocess,change(i),degconst*angledelta(i));
                        else
                            fprintf(...
                                'Dataset %d, step_IVA %d - IVA_lambda %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\n',...
                                i,step_IVA,lambda_IVA_preprocess,change(i),degconst*angledelta(i),(ncomps(i)-sum(diag(signs{i})))/2);
                        end
                    elseif ~extended
                        fprintf(...
                            'Dataset %d, step_IVA %d - IVA_lambda %5f, wchange %7.6f\n',i,step_IVA,lrate(i),change(i));
                    else
                        fprintf(...
                            'Dataset %d, step_IVA %d - IVA_lambda %5f, wchange %7.6f, %d subgauss\n',...
                            i,step_IVA,lambda_IVA_preprocess,change(i),(ncomps(i)-sum(diag(signs{i})))/2);
                    end % step_IVA > 2
                end; % if verbose
            end
            
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate if degree is larger than 60%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if  (~isempty(find(degconst*angledelta > annealdeg))) %%%(degconst*angledelta(1) > annealdeg) |(degconst*angledelta(2) > annealdeg)|(degconst*angledelta(3) > annealdeg) %%degconst*angledelta_all > annealdeg %%entropychange(i)<0 %|
                olddelta  = delta;                % accumulate angledelta until
                oldchange  = change;              %  annealdeg is reached
                olddelta_concat = delta_concat;
                oldchange_all = change_all;
                lambda_IVA_preprocess = lambda_IVA_preprocess*annealstep; % anneal learning rate                
            elseif step_IVA == 1                     % on first step_IVA only
                olddelta   = delta;                  % initialize
                oldchange  = change;
                olddelta_concat = delta_concat;
                oldchange_all = change_all;
            end
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            if (step_IVA >2) & (length(find(change < nochange)) == length(data))
                for i = 1:length(data)   
                    stopsign_IVA(i)=1;   % stop when weights stabilize
                end
            elseif step_IVA >= maxsteps
                for i = 1:length(data)   
                    stopsign_IVA(i)=1;   % stop when it reaches the maximum steps 
                end                
            elseif (~isempty(find(change > DEFAULT_BLOWUP)))%%%(change(1) > DEFAULT_BLOWUP)|(change(2) > DEFAULT_BLOWUP)|(change(3) > DEFAULT_BLOWUP) % if weights blow up,
                lambda_IVA_preprocess = lambda_IVA_preprocess*DEFAULT_BLOWUP_FAC; % with lower learning rate
            end  
            %%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            oldweights = weights;           
        end % end if weights in bounds        
    end
end



              
         