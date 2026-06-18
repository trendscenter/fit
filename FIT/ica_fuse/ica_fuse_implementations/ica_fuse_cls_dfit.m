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

