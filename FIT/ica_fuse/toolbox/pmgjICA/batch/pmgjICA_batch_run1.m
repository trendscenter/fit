numofSub=100;
dim_modality=[53, 63, 52];
num_fmri=53;
preProceType=5;
num_mod=2;
num_pc1=16; 
num_pc2=10; 

save_results='/data/users4/ceierud/proj25/misc/ibrahim022625/pmlJICA_fromgift/res18sub100';
save_model_info='/data/users4/ceierud/proj25/misc/ibrahim022625/pmlJICA_fromgift/res18sub100/tst_ica_parameter_info.mat';

data_path='/data/users2/ibrahim/pmlJICA_fromgift/preprocessed_fbirn_notlc_rmpt_sep.mat';
gm_mask_path='/data/users4/ceierud/proj25/misc/ibrahim022625/code/smre_gm_mark042425sub100.mat';
fmri_mask_path = '/data/users4/ceierud/proj25/misc/ibrahim022625/code/icn_mask042425sub100.mat';
