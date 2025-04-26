
function [weights,sphereGmCm, datacombined_pca] = mygift_pmljICA(dataMat, sesInfo, ICA_Options)
    
    if isfield(sesInfo, 'outputDir')
        outputDir = sesInfo.outputDir;
    else
        outputDir=sesInfo.userInput.pwd;
    end
    pcaType=sesInfo.pcaType;

    %%%size(whitesigGm) %% 74994          100
    %%size(dewhiteMGm)  %%% 13780          100
    
    %%size(dataMat)  %%% 15      143733
  
    if (sesInfo.userInput.groupica_algorithm==1)
        disp('Inside new algo.........')
        if sesInfo.userInput.ICAcallcount==1
            disp('First.........')
            numOfPC=sesInfo.reduction(1).numOfPCAfterReduction;
            numOfPC
        else
            numOfPC=sesInfo.reduction(2).numOfPCAfterReduction;
            numOfPC
        end
        pcaType=sesInfo.pcaType;
        dataGm1=dataMat(:, 1:length(sesInfo.userInput.maskmod1)); %%% subjectxVoxels      
        dataCm1=dataMat(:, length(sesInfo.userInput.maskmod1)+1:end);  %%%whitesigCm';
        %% Do PCA with whitening:  Input: VoxelXSubejcts
        [pcasig, dewhiteM, Lambda, V, whiteM] = mygift_calculate_pca(dataGm1', numOfPC, 'type', pcaType); 
        dataGm=pcasig';
        if sesInfo.userInput.ICAcallcount==2
            pcaout=fullfile(outputDir, [sesInfo.data_reduction_mat_file, num2str(sesInfo.userInput.ICAcallcount), '-1_', num2str(1),'.mat']);
            mygift_save(pcaout, 'pcasig', 'dewhiteM', 'Lambda', 'V', 'whiteM');
        end
        [pcasig2, dewhiteM2, Lambda2, V2, whiteM2] = mygift_calculate_pca(dataCm1', numOfPC, 'type', pcaType); 
        dataCm=pcasig2'; 

        if sesInfo.userInput.ICAcallcount==1
            pcasig_sepcom=[pcasig; pcasig2];
            deWhiteM_sepcom=[dewhiteM;dewhiteM2];
            whiteM_sepcom=[whiteM;whiteM2];
            pcaout=[sesInfo.data_reduction_mat_file, '_sep_combined', '.mat'];           
       
            pcajoint = fullfile(outputDir, pcaout);                              
            mygift_save(pcajoint, 'pcasig_sepcom', 'deWhiteM_sepcom', 'whiteM_sepcom');
        end


        if sesInfo.userInput.ICAcallcount==2
            pcasig=pcasig2;
            dewhiteM=dewhiteM2;
            Lambda=Lambda2;
            V=V2;
            whiteM=whiteM2;
            pcaout=fullfile(outputDir, [sesInfo.data_reduction_mat_file, num2str(sesInfo.userInput.ICAcallcount), '-1_', num2str(2),'.mat']);
            mygift_save(pcaout, 'pcasig', 'dewhiteM', 'Lambda', 'V', 'whiteM');
        end 

        
        clear pcasig dewhiteM Lambda V whiteM pcasig2 dewhiteM2 Lambda2 V2 whiteM2;

    else
        numOfPC=sesInfo.reduction(2).numOfPCAfterReduction;
        %%% dataMat: nCompsxVoxels
        dataGm=dataMat(:, 1:length(sesInfo.userInput.maskmod1));       
        dataCm=dataMat(:, length(sesInfo.userInput.maskmod1)+1:end);  %%%whitesigCm';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('check data..')
    size(dataGm)
    size(dataCm)

    [chansCm frames] = size(dataCm); % determine the data size
    %DEFAULT_LRATE = 0.015/log(chans); 
    lrate = 0.015/log(chansCm);
    
    %%%%%%%% Initialize weights 
    if (sesInfo.userInput.groupica_algorithm==1)
        ncomps=chansCm;
    else
        ncomps=sesInfo.ICA_Options{20};
    end
    weights = eye(ncomps,chansCm); %%% begin with the identity matrix Initilize here first    
    %%%weights = randsmall(ncomps,chansGm);
    size(weights)
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_mods=2;
    pass=1;
    max_steps=500;  %%% 2000

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% imortant for termination and optimization %%%%%%%%
    degconst = 180./pi;
    annealstep = 0.9;%%
    %%%annealdeg  = DEFAULT_ANNEALDEG - momentum*90; % heuristic
    annealdeg  = 60;
    %%%%%%% imortant for termination and optimization %%%%%%%%
    
    
    nochange=0.00004; %0.0000001;
    DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
    DEFAULT_BLOWUP_FAC   = 0.9;
    posactflag='off';
    wts_passed=0;
    
    %%%load weightsdiff_700.mat;
    
    fg=figure;
    
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
        %CE042425 error here
        [W1, lrate, ~] = pmljICA_optimize(dataGm, maxsteps, lrate, weights, ICA_Options{1:length(ICA_Options)});  %%%% verbose=0;  %%% always zero, need change inside if...
        weights=W1; 
        
        %%disp('weights...')
        %%size(weights)  %%% 15 15
    
        [W2, lrate, sphere, data, signs, bias] = pmljICA_optimize(dataCm, maxsteps, lrate, weights, ICA_Options{1:length(ICA_Options)});
        %%weights=W2;
        
        %%size(weights)
        %%%%oldwtchange=W2-W1;
        weights=(W1+W2)./num_mods;
        
        %disp('learning rate before update in the outer loop....')
        %lrate
    
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
                   annealstep = 0.95;%%               
               end
               
               lrate = lrate*annealstep;          % anneal learning rate
               olddelta   = delta;                % accumulate angledelta until
               oldchange  = change;               %  annealdeg is reached                      
            end
    
            fprintf(...
                      'step %d - lrate %7.9f, wchange %7.7f, angledelta %4.1f deg\n',...
                      pass,lrate,change,degconst*angledelta);
    
            
    
            oldweights=weights;
            
            if pass==2
                plt_change(pass)=change;
                plot(2:pass, plt_change(2:pass), '-ob');
                %xlim([1,max_steps+1]);
                %xticks(1:1:max_steps+1);     
                set(gca,'XTick',(1:50:701))
                %%hold on
                %%pause(0.01)
                grid on
                title('Jointly weights optimization')
                xlabel('Epoch')
                ylabel('Weights')
                hold on
            elseif pass<max_steps
                %%%%plt_change(step)=change;
                plt_change(pass)=change;
                plot(2:pass, plt_change(2:pass), '-ob');
                %xlim([1,max_steps+1]);
                %xticks(1:1:max_steps+1);     
                set(gca,'XTick',(1:50:701))
                %%hold on
                pause(0.01)
                grid on
                title('Jointly Optimization')
                xlabel('Epoch')
                ylabel('Weight change')
            end
    
        end
        
        pass=pass+1;    
        %
        %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %disp('pass and change')
        %pass
        %change
        
        if pass >2 & change < nochange,      % apply stopping rule
           laststep=pass;            
           pass=max_steps+1;                  % stop when weights stabilize
        elseif change > DEFAULT_BLOWUP,      % if weights blow up,
           lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying 
        end; 
    
    
    end
    hold off;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
    
     %%%%pcasig=whitesigjoint
    if (sesInfo.userInput.groupica_algorithm==1 || strcmpi(sesInfo.userInput.pjICA, 'yes'))

        [pcasigJoint, dewhiteMJoint, LambdaJoint, VJoint, whiteMJoint] = mygift_calculate_pca(dataMat', numOfPC, 'type', pcaType); 
        sphereGmCm = 2.0*inv(sqrtm(cov(pcasigJoint))); % find the "sphering" matrix = spher()  %%% No need sphereing, since we have used in PCA
        sphereGm = 2.0*inv(sqrtm(cov(dataGm')));  %%% voxelsxcomp
        sphereCm = 2.0*inv(sqrtm(cov(dataCm')));

        %%%% whiteningMatrix = sqrtm(Lambda) \ V';
        %%%% dewhiteningMatrix = V * sqrtm(Lambda);

        %%size(dataGm')  %%% voxelsxcomp

        %%%% sphereGmCm: 30x30
        size(sphereGmCm)
        %%datacombined_pca=pcasigJoint';  %% dataGmCm
        datacombined_pca=[dataGm, dataCm];
        %%%%%%%%%%%%%% Orient components towards positive activation %%%%%%%%%%%
        dataj = sphereGmCm*(pcasigJoint');
    
    else
        jointFile=fullfile(sesInfo.userInput.pwd, "joint_pcaofICAdata.mat");
        load(jointFile);    
        sphereGmCm = 2.0*inv(sqrtm(cov(whitesigjoint))); % find the "sphering" matrix = spher()  %%% No need sphereing, since we have used in PCA
        %%%% sphereGmCm: 30x30
        size(sphereGmCm)
        %%datacombined_pca=whitesigjoint';  %%%% dataGmCm
        datacombined_pca=[dataGm, dataCm];
        %%%%%%%%%%%%%% Orient components towards positive activation %%%%%%%%%%%
        dataj = sphereGmCm*(whitesigjoint');
    end
    size(datacombined_pca)
    disp('dataj...')
    size(dataj)

    if strcmp(posactflag,'on')
       [activations, winvout, weights] = mygift_posact(dataj, weights);
       % changes signs of activations and weights to make activations
       % net rms-positive
    else
       activations = weights*dataj;
    end
    
    fprintf(...
          'Sorting components in descending order of mean projected variance ...\n');    
    
    if wts_passed == 0
       %
       %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %
       meanvar  = zeros(ncomps,1);      % size of the projections
       if ncomps == chansCm % if weights are square . . . %% equal to number of component
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
       %%bias  = bias(windex);		% reorder them
       %%%signs = diag(signs);        % vectorize the signs matrix
       %%%signs = signs(windex);      % reorder them
    else
       fprintf('Components not ordered by variance.\n');
    end    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (sesInfo.userInput.groupica_algorithm==1)
        if sesInfo.userInput.ICAcallcount==2
            pcaout=[sesInfo.data_reduction_mat_file, num2str(sesInfo.userInput.ICAcallcount), '-1_', 'joint','.mat'];
        else 
            pcaout=[sesInfo.data_reduction_mat_file, '_firststage_ref', '.mat'];           
        end
        pcajoint = fullfile(outputDir, pcaout);                              
        mygift_save(pcajoint, 'pcasigJoint', 'dewhiteMJoint', 'LambdaJoint', 'VJoint', 'whiteMJoint'); 
        close(fg);
        %%clf(fg);
    elseif strcmpi(sesInfo.userInput.pjICA, 'yes')
        ica_data_path = fullfile(outputDir, 'sphere_weights.mat'); 
        mygift_save(ica_data_path, 'sphereGmCm', 'sphereGm', 'sphereCm', 'weights');
    end
    
%     W1=weights*sphereGmCm;   %%%% 30x30
%     sources=W1*dataGmCm;
%     A = pinv(W1);  %%%% A30sp: 30x30 : should be used 
%     
%     
%     A = dewhiteMGmCm*A;  %%%% dewhiteMGmCm: 13780x30
%     W = pinv(A);

