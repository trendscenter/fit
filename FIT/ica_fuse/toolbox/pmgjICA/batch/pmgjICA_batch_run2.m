num_subjects=100;
hc_sub=50;
ad_sub=50;
num_fmri=53;
sfilename='';


source='/data/users4/ceierud/proj25/misc/ibrahim022625/pmlJICA_fromgift/res18sub100/source042525';

procedure='dual_regress_zmean_joint_notzscoredat2'; %%1

param_pmlg = '/data/users4/ceierud/proj25/misc/ibrahim022625/pmlJICA_fromgift/res18sub100/tst_ica_parameter_info.mat';


convert_to_zscores='yes'; %%% it is needed
thr=1.5; %ce042525 was 0.01;  %%%% first one threshold and min/max range
% i_ad=50;

p_val=[];





