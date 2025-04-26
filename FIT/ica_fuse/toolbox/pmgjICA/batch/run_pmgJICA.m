function run_pmgjICA(s_batch)

%%Parallel Multilink Group Joint ICA

%%%% Parameter setup
stru_in = ica_fuse_eval_script(s_batch);

pmgJICA_recon(stru_in.numofSub, stru_in.save_model_info, stru_in.num_pc1, stru_in.num_pc2, stru_in.num_mod, stru_in.num_fmri, stru_in.dim_modality, stru_in.preProceType,...
    stru_in.data_path, stru_in.save_results, stru_in.gm_mask_path, stru_in.fmri_mask_path)

disp('Done with 1st step of Parallel Multilink Group Joint ICA')
%%  please see the output in the "save_results"
