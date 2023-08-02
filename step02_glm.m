clear
clc

load e.mat

model_name = 'model_01';

dirStats = e.mkdir('glm',model_name);
dirFunc  = get_parent_path( e.getSerie('fmri_nm_bold_Run01').getVolume('s5wts').toJob() );
dirFunc   = cellfun(@cellstr, dirFunc, 'UniformOutput', 0);

onsetspath = '/home/benoit.beranger/tmp/MOTIV_APP/onsets';
e.getSerie('fmri_nm_bold_Run01').addStim(onsetspath, 'onsets.mat', model_name)
onsets = e.getSerie('fmri_nm_bold').getStim(model_name).toJob(0);

clear par
par.file_reg = '^s5wts_.*nii';
par.rp       = 1;
par.rp_regex = '^multiple_regressors.txt';

% Masking
par.mask_thr = 0.1; % spm default option
par.mask     =  {}; % cell(char) of the path for the mask of EACH model : N models means N paths

par.sge      = 0;
par.run      = 1;
par.display  = 0;
par.redo     = 0;

par.TR = 1.660;

job_first_level_specify(dirFunc,dirStats,onsets,par);
e.addModel('glm',model_name,model_name);
mdl = e.getModel(model_name);
fspm = mdl.getPath();


%%

clear par
par.write_residuals = 0;

par.jobname  = 'spm_glm_est';
par.walltime = '11:00:00';

par.sge      = 0;
par.run      = 1;
par.display  = 0;
par.redo     = 0;

job_first_level_estimate(fspm, par);


%%

clear par
par.sessrep         = 'none';
par.report          = 0;

par.jobname         ='spm_glm_con';
par.walltime        = '04:00:00';

par.sge             = 0;
par.run             = 1;
par.display         = 0;
par.delete_previous = 1;


learn  = [1 0];
recall = [0 1];


contrast_T.values = {

learn
recall

learn - recall
recall- learn

}';

contrast_T.names = {

'learn'
'recall'

'learn - recall'
'recall- learn'

}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));

contrast_F.names  = {'F-all'};
contrast_F.values = {eye(2)};
contrast_F.types  = cat(1,repmat({'F'},[1 length(contrast_F.names)]));

contrast.names  = [contrast_F.names  contrast_T.names ];
contrast.values = [contrast_F.values contrast_T.values];
contrast.types  = [contrast_F.types  contrast_T.types ];

job_first_level_contrast(fspm,contrast,par);


%%

mdl.show()
