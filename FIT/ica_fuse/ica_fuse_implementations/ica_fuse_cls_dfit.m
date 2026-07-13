classdef ica_fuse_cls_dfit
    % Cyrus Eierud Date 6/16/26
    % Code that supports FIT processing in FIT from GIFT dynamic results
    % (e.g., dynamic FNC)
    % Example to instantiate dfit:
    %       oc_dfit = ica_fuse_cls_dfit('/path/to/gift_param.mat');
    
    properties
        clusterInfo
        s_fit_outputdir
        s_prefix
        s_gift_path
        b_dfit_selected_in_batch_or_gui
    end
    
    methods       
        function n_ret = dyn_states_gift2fit_save(obj)
            % Converts GIFT data for FIT usage;
            n_ret = 1;
            ix_subj = 1:size(obj.clusterInfo.states,1);
            ix_sess = 1:size(obj.clusterInfo.states,2);
            ix_states = 1:size(obj.clusterInfo.Call,1);
            for n_subj = ix_subj
                for n_sess = ix_sess
                    tmp = load([obj.s_gift_path filesep obj.s_prefix ...
                        '_dfnc_sub_' sprintf('%03d', n_subj) '_sess_' ...
                        sprintf('%03d', n_sess) '_results']);
                    con_states = squeeze(obj.clusterInfo.states(n_subj,n_sess,:));
                    for n_state = ix_states
                        ix_state = find(con_states == n_state);
                        if isempty(ix_state)
                            %ce061626 count these and disclain how many
                            %approximations were done
                            % Subject not in this state, substitute with centroid
                            FNCdyn = obj.clusterInfo.Call(n_state,:)';
                        else
                            FNCdyn = mean(tmp.FNCdyn(ix_state,:))';
                        end
                        s_file_save = [obj.s_fit_outputdir filesep ...
                            'dfit' filesep obj.s_prefix ...
                            '_sub-' sprintf('%03d', n_subj) ...
                            '_sess-' sprintf('%03d', n_sess) ...
                            '_kmnstate-' sprintf('%02d', n_state) '_avgConn.mat'];
                        ica_fuse_save(s_file_save, 'FNCdyn'); 
                        clear FNCdyn ix_state;
                    end
                end
            end

            n_ret = 0; %0=worked
        end
        
        function obj = ica_fuse_cls_dfit(s_input_file, s_fit_outputdir)
            % Constructor, load valid parameter file

            obj.b_dfit_selected_in_batch_or_gui = 1;

            try
                stru_tmp = ica_fuse_read_variables(s_input_file, 'dynfitFileGiftDfncPostProcessMat', {'character'});
            catch
                % Assuming that dynfitFileGiftDfncPostProcessMat did not
                % exist in {prefix}_dfnc_post_process.mat, so skip
                % everything dynamic
                obj.b_dfit_selected_in_batch_or_gui = 0;
                return
            end

            %Simple validation and initiation
            if exist(stru_tmp.dynfitFileGiftDfncPostProcessMat,'file')
                drawnow; %File found, just continue
            else
                error(['GIFT {prefix}_dfnc_post_process.mat file for dynamic FIT not found:' s_gift_param]);
            end           
            try
                load(stru_tmp.dynfitFileGiftDfncPostProcessMat);
                obj.clusterInfo = clusterInfo;
                clear clusterInfo FNCamp FNCcm meta_states_info sgica;
            catch
                error(['Problems reading GIFT file for dynamic FIT:' s_gift_param]);
            end

            [s_gift_path, s_name, s_ext] = fileparts(stru_tmp.dynfitFileGiftDfncPostProcessMat);
            s_suffix = lower('_dfnc_post_process');
            if endsWith(lower(s_name), s_suffix)
                s_prefix = s_name(1:end-length(s_suffix));
                obj.s_prefix = s_prefix;
                obj.s_gift_path = s_gift_path;
            else
                error('dynfitFileGiftDfncPostProcessMat does not end with: _dfnc_post_process.mat');
            end


            if exist(s_fit_outputdir,'dir')
                % Create destination folder for dfit data
                try
                    mkdir([s_fit_outputdir filesep 'dfit']);
                catch
                end
                obj.s_fit_outputdir = s_fit_outputdir;
            else
                error(['FIT output folder not found:' s_fit_outputdir]);
            end               
        end

        function n_ret = dyn_matching_components_across_states(obj)
            n_states = size(obj.clusterInfo.Call,1);
            
            tmp_files = dir(fullfile(obj.s_fit_outputdir, [obj.s_prefix '_state1_*comp_br_comb_1.mat']));                                 
            s_tmp = load([tmp_files(1).folder filesep tmp_files(1).name]);
% %             clear tmp_files
% %             [obj.s_prefix  '_state' num2str(1) '_*_' num2str(n_comps) 'comp_joint_comp_ica_feature_2_' sprintf('%03d', comp1) '.asc']
            
            n_comps = size(s_tmp.icasig,1);
            clear s_tmp tmp_files;
            match_modality = "struct"; % options: "struct", "dFNC"
            outpath = [obj.s_fit_outputdir filesep 'post_processing/'];
            try
                mkdir outpath;
            catch
            end

            % all pairs of states to perform cross-state component matching
            state_pairs = nchoosek(1:n_states, 2);

            % comp_matches - size [n_state_pairs x n_comps x 2]
            % contains idxs for best matched component pairs
            % for each state pair in state_pairs (above)
            comp_matches = zeros(size(state_pairs,1),n_comps,2);

            % matched_corrs - size [n_pairs]
            matched_corrs = zeros(size(state_pairs,1),n_comps);

            % option 1 - find cross-fusion component matches based on structural components
            % most common - we want to find the best matched components across state fusions
            % to identify "static" or "dynamic" structural patterns

                 for pair = 1:size(state_pairs,1)

                    state1 = state_pairs(pair,1); state2 = state_pairs(pair,2);

                    corrs = zeros(n_comps,n_comps);

                    % loop through all comps in state 1
                    for comp1 = 1:n_comps

                        % feature_2 is the structural feature based on our batch file implementation
                        clear tmp_files; s_file_wild_card=[obj.s_prefix  '_state' num2str(state1) '_*_' num2str(n_comps) 'comp_joint_comp_ica_feature_2_' sprintf('%03d', comp1) '.asc'];
                        tmp_files = dir(fullfile(obj.s_fit_outputdir, s_file_wild_card));                                                       
                        if ~(size(tmp_files,1) == 1)
                            error(['Error in ica_fuse_cls_dfit: exactly one file should match ' obj.s_fit_outputdir filesep s_file_wild_card]);
                        end                        
%                         s1 = load(append("demo_results/dynamicFusion_demo_dFNC_state",string(state1),"_GMV_5comp_joint_comp_ica_feature_2_", sprintf('%03d', comp1),".asc"));
                        s1 = load([tmp_files(1).folder filesep tmp_files(1).name]);
                        % loop through all comps in state 2
                        for comp2 = 1:n_comps
                            clear tmp_files; s_file_wild_card=[obj.s_prefix  '_state' num2str(state2) '_*_' num2str(n_comps) 'comp_joint_comp_ica_feature_2_' sprintf('%03d', comp2) '.asc'];
                            tmp_files = dir(fullfile(obj.s_fit_outputdir, s_file_wild_card));                                            
                            if ~(size(tmp_files,1) == 1)
                                error(['Error in ica_fuse_cls_dfit: exactly one file should match ' obj.s_fit_outputdir filesep s_file_wild_card]);
                            end                        
%                             s2 = load(append("demo_results/dynamicFusion_demo_dFNC_state",string(state2),"_GMV_5comp_joint_comp_ica_feature_2_", sprintf('%03d', comp2),".asc"));
                            s2 = load([tmp_files(1).folder filesep tmp_files(1).name]);

                            disp(append("Computing state matches: State ", string(state1)," Comp ", string(comp1), " & State ", string(state2), " Comp ", string(comp2)));

                            % compute the correlation
                            c = corrcoef(s1(:,2), s2(:,2));

                            % fill in correlation in corrs matrix
                            corrs(comp1,comp2) = c(1,2);

                        end
                    end

                    % initializations
                    idxs = zeros(n_comps,2);
                    c = abs(corrs);
                    i = 1;

                    % greedy matching scheme - find highest abs(corr) pair of components
                    % record pair in idxs
                    % zero out that row/col to remove from consideration at next iteration
                    while sum(c,"all") > 0

                        [row, col] = find(c == max(c, [], "all"));
                        idxs(i,:) = [row, col];
                        c(row,:) = zeros(1,n_comps); c(:,col) = zeros(n_comps,1);
                        i = i+1;

                    end

                    % sort corrs by idxs of our best matched pairs
                    corrs = corrs(idxs(:,1),idxs(:,2));

                    % the diagonal of this sorted matrix is now the correlations
                    % for the best matched pairs
                    disp([idxs, diag(corrs)]);

                    % comp_matches records the matched component pairs 
                    comp_matches(pair, :, :) = idxs;

                    % matched_corrs records the correlations for those matched pairs
                    matched_corrs(pair,:) = diag(corrs);

                end

                % save results
                save(append(outpath,"dynamicFusion_postprocessing_component_matches_struct.mat"), "comp_matches", "matched_corrs");

            % option 2 - find cross-fusion component matches based on dFNC
            % not common - we don't expect much commonality between the dFNC components
            % across state fusions because by definition these are different FNC states


                 for pair = 1:size(state_pairs,1)

                    state1 = state_pairs(pair,1); state2 = state_pairs(pair,2);

                    corrs = zeros(n_comps,n_comps);

                    % loop through all comps in state 1
                    for comp1 = 1:n_comps

                        % feature_2 is the structural feature based on our batch file implementation
                        clear tmp_files; s_file_wild_card=[obj.s_prefix  '_state' num2str(state1) '_*_' num2str(n_comps) 'comp_joint_comp_ica_feature_1_' sprintf('%03d', comp1) '.asc'];
                        tmp_files = dir(fullfile(obj.s_fit_outputdir, s_file_wild_card));                                                       
                        if ~(size(tmp_files,1) == 1)
                            error(['Error in ica_fuse_cls_dfit: exactly one file should match ' obj.s_fit_outputdir filesep s_file_wild_card]);
                        end                        
                        s1 = load([tmp_files(1).folder filesep tmp_files(1).name]);
                        % loop through all comps in state 2
                        for comp2 = 1:n_comps
                            clear tmp_files; s_file_wild_card=[obj.s_prefix  '_state' num2str(state2) '_*_' num2str(n_comps) 'comp_joint_comp_ica_feature_1_' sprintf('%03d', comp2) '.asc'];
                            tmp_files = dir(fullfile(obj.s_fit_outputdir, s_file_wild_card));                                            
                            if ~(size(tmp_files,1) == 1)
                                error(['Error in ica_fuse_cls_dfit: exactly one file should match ' obj.s_fit_outputdir filesep s_file_wild_card]);
                            end                        
                            s2 = load([tmp_files(1).folder filesep tmp_files(1).name]);

                            disp(append("Computing state matches: State ", string(state1)," Comp ", string(comp1), " & State ", string(state2), " Comp ", string(comp2)));

                            % compute the correlation
                            c = corrcoef(s1(:,2), s2(:,2));

                            % fill in correlation in corrs matrix
                            corrs(comp1,comp2) = c(1,2);

                        end
                    end

            
            
            
            
            

                    % initializations
                    idxs = zeros(n_comps,2);
                    c = abs(corrs);
                    i = 1;

                    % greedy matching scheme - find highest abs(corr) pair of components
                    % record pair in idxs
                    % zero out that row/col to remove from consideration at next iteration
                    while sum(c,"all") > 0

                        [row, col] = find(c == max(c, [], "all"));
                        idxs(i,:) = [row, col];
                        c(row,:) = zeros(1,n_comps); c(:,col) = zeros(n_comps,1);
                        i = i+1;

                    end

                    % sort corrs by idxs of our best matched pairs
                    corrs = corrs(idxs(:,1),idxs(:,2));

                    % the diagonal of this sorted matrix is now the correlations
                    % for the best matched pairs
                    disp([idxs, diag(corrs)]);

                    % comp_matches records the matched component pairs 
                    comp_matches(pair, :, :) = idxs;

                    % matched_corrs records the correlations for those matched pairs
                    matched_corrs(pair,:) = diag(corrs);

                end

                % save results
                save(append(outpath,"dynamicFusion_postprocessing_component_matches_dFNC.mat"), "comp_matches", "matched_corrs");
        end        
        
        function obj = set_b_dfit_selected_in_batch_or_gui(obj, newValue)
            % Method to set Value
            obj.b_dfit_selected_in_batch_or_gui = newValue;
        end

        function val = get_b_dfit_selected_in_batch_or_gui(obj)
            % Method to get Value
            val = obj.b_dfit_selected_in_batch_or_gui;
        end        
    end
end

