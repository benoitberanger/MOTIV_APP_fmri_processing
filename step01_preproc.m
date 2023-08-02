%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   MULTI - ECHO    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear


%% Requirements

% what : matvol <==== CENIR toolbox for MRI analysis, it contains all wrappers for all libraries (AFNI, FSL, ANTs,....)
% language : MATLAB
% where : https://github.com/romainVala/matvol
% comments : must be in MATLAB path
assert( ~isempty(which('init_path_matvol')), 'matvol not found' )

% what : SPM12
% language : MATLAB
% where : https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
% comments : must be in MATLAB path
% maybe need to initialize SPM paths with : spm_jobman('initcfg')
assert( ~isempty(which('spm_read_vols')), 'SPM12 not found' )

% what : CAT12
% language : MATLAB + SPM12
% where : https://neuro-jena.github.io/cat/
% comments : must be in spm/toolbox/ dir
% maybe need to initialize SPM paths with : spm_jobman('initcfg')
assert( ~isempty(which('cat12')), 'cat12 not found' )

% what : TAPAS/PhysIO
% language : MATLAB + SPM12
% where : https://github.com/translationalneuromodeling/tapas/tree/master/PhysIO
% comments : must be in spm/toolbox/ dir
% maybe need to initialize SPM paths with : spm_jobman('initcfg')
assert( ~isempty(which('tapas_physio_init')), 'tapas_physio_init not found' )

% what : CMRR_MB
% language : MATLAB
% where : https://github.com/CMRR-C2P/MB
% comments : must be in MATLAB path
assert( ~isempty(which('readCMRRPhysio')), 'readCMRRPhysio not found' )

% what : AFNI
% language : C
% where : https://afni.nimh.nih.gov/
% comments : must be detectable my MATLAB underlysing terminal
assert( unix('which afni_proc.py') == 0, 'afni not found' )

% what : FSL
% language : C
% where : https://fsl.fmrib.ox.ac.uk/fsl/fslwiki
% comments : must be detectable my MATLAB underlysing terminal
assert( unix('which fslmaths') == 0, 'fslmaths not found' )


%% main variables

main_dir  = '/home/benoit.beranger/tmp/MOTIV_APP';
nifti_dir = fullfile(main_dir,'nifti');

% local or cluster ?
global_parameters.run = 1;
global_parameters.sge = 0;


%% Build exam lists

% get exam==subject (directory) list
e = exam(nifti_dir, 'MOTIV_APP_Pilote\d{2}'); % all subjects with multi-echo

% add T1 serie (directory) in each exam
e.addSerie('T1w$', 'anat_T1', 1 );
e.getSerie('anat').addVolume('^v_.*nii', 'v', 1);

run_list = {'Run01', 'Run02'};
for r = 1 : length(run_list)
    run_name = run_list{r};
    e.addSerie([run_name '$'            ] , ['fmri_nm_bold_' run_name], 1 );
    e.addSerie([run_name '_revPE$'      ] , ['fmri_rv_bold_' run_name], 1 );
    e.addSerie([run_name '_SBRef$'      ] , ['fmri_nm_sbref_' run_name], 1 );
    e.addSerie([run_name '_revPE_SBRef$'] , ['fmri_rv_sbref_' run_name], 1 );
    e.addSerie([run_name '_PhysioLog$'  ] , ['phy_nm_' run_name], 1 );
end

e.getSerie('fmri_.*_bold').addVolume('^v_.*nii$'    , 'v', 3);
e.getSerie('fmri_.*_sbref').addVolume('^v_.*e1.nii$', 'v', 1);
e.getSerie('phy').addPhysio('dcm$', 'dcm', 1);


%% check if all echos have the same number of volumes #MATLAB/matvol
% this step takes time the first time you run it, but after its neglectable

e.getSerie('fmri').getVolume('^v').removeEmpty().check_multiecho_Nvol()


%% segement T1 #MATLAB/SPM::CAT12

anat = e.gser('anat_T1').gvol('^v');
job_do_segmentCAT12(anat,global_parameters);


%% Sort echos #MATLAB/matvol

cfg = struct;
cfg.run  = 1;
cfg.fake = 0;
cfg.sge  = 0;
cfg.redo = 0;

meinfo = job_sort_echos( e.getSerie('fmri_nm_bold') , cfg );


%% minimal preprocessing for multi-echo #bash/ANFI::afni_proc.py

cfg = global_parameters;
cfg.blocks  = {'tshift', 'volreg', 'blip'};
afni_prefix = char(cfg.blocks); % {'tshift', 'volreg', 'blip'}
afni_prefix = afni_prefix(:,1)';
afni_prefix = fliplr(afni_prefix); % 'bvt'
afni_subdir = ['afni_' afni_prefix];
cfg.subdir = afni_subdir;

cfg.blip.forward = e.getSerie('fmri_nm_sbref').getVolume('^v');
cfg.blip.reverse = e.getSerie('fmri_rv_sbref').getVolume('^v');

job_afni_proc_multi_echo( meinfo, cfg );


%% mask EPI #bash/FSL

fin  = e.getSerie('fmri_nm_bold').getVolume(['^' afni_prefix 'e1']);
do_fsl_robust_mask_epi( fin, global_parameters );

% Checkpoint & unzip
cfg = global_parameters;
cfg.jobname = 'unzip_and_keep__bet';
e.getSerie('fmri_nm_bold').getVolume(['^bet_Tmean_' afni_prefix 'e1$']).removeEmpty().unzip_and_keep(cfg)


%% echo combination #python/TEDANA

tedana_subdir = ['tedana0011_' afni_prefix];
job_tedana_0011( meinfo, afni_prefix, tedana_subdir, ['bet_Tmean_' afni_prefix 'e1_mask.nii.gz'], global_parameters );

% Checkpoint & unzip
cfg = global_parameters;
cfg.jobname = 'unzip_and_keep__tedana';
e.getSerie('fmri_nm_bold').getVolume('^ts_OC').removeEmpty().unzip_and_keep(cfg)


%% coregister EPI to anat #MATLAB/SPM12

cfg = global_parameters;
cfg.type  = 'estimate';

src = e.getSerie('fmri_nm_bold').removeEmpty().getVolume(['^bet_Tmean_' afni_prefix 'e1$']);
oth = e.getSerie('fmri_nm_bold').removeEmpty().getVolume('^ts_OC');
ref = e.getSerie('fmri_nm_bold').removeEmpty().getExam.getSerie('anat_T1').getVolume('^p0');

cfg.jobname = 'spm_coreg_epi2anat';
job_coregister(src,ref,oth,cfg);


%% normalize EPI to MNI space #MATLAB/SPM12

img4D = e.getSerie('fmri_nm_bold').getVolume('^ts_OC').removeEmpty();
img3D = e.getSerie('fmri_nm_bold').getVolume(['^bet_Tmean_' afni_prefix 'e1$']).removeEmpty();

img = img4D + img3D;
y   = img.getExam.getSerie('anat_T1').getVolume('^y');
cfg = global_parameters;
if global_parameters.sge % for cluster all jobs preparation, we need to give the voxel size
    % !!! this assumes all runs have the same resolution !!!
    % fetch resolution from the first volume
    V = e.getSerie('fmri_nm_bold').getVolume('^v');
    V = spm_vol(deblank(V(1).path(1,:)));
    cfg.vox = sqrt(sum(V(1).mat(1:3,1:3).^2));
end
job_apply_normalize(y,img,cfg);


%% smooth EPI #MATLAB/SPM12

cfg = global_parameters;
img = e.getSerie('fmri_nm_bold').getVolume('^w.*_OC').removeEmpty();

cfg.smooth   = [5 5 5];
cfg.prefix   = 's5';
cfg.jobname  = 'spm_smooth5';
job_smooth(img,cfg);


%% coregister WM & CSF on functionnal (using the warped mean) #SPM12
% This will be used for TAPAS:PhysIO

cfg = global_parameters;
ref = e.getSerie('fmri_nm_bold');
ref = ref(:,1).getVolume(['wbet_Tmean_' afni_prefix 'e1']);
src = e.getSerie('anat_T1').getVolume('^wp2');
oth = e.getSerie('anat_T1').getVolume('^wp3');
cfg.type = 'estimate_and_write';
cfg.jobname = 'spm_coreg_WMCSF2wEPI';
job_coregister(src,ref,oth,cfg);


%% rp afni2spm #matlab/matvol

% input
dfile = e.getSerie('fmri_nm_bold').getRP('rp_afni').removeEmpty();

% output
output_dir = fullfile( dfile.getSerie().getPath(), tedana_subdir );

% go
job_rp_afni2spm(dfile, output_dir);


%% extract physio from special dicom

% https://github.com/CMRR-C2P/MB

e.getSerie('phy').getPhysio('dcm').extract()

% e.getSerie('phy').getPhysio('phy').check() % takes a bit of time, use it once to verify your data


%% PhysIO nuisance regressor generation #matlab/TAPAS-PhysIO
%% Prepare files

% get physio files & check if some are missing
info = e.getSerie('phy').removeEmpty().getPhysio('info');   info = info(:);   missing_info = cellfun( 'isempty', info(:).getPath() );
puls = e.getSerie('phy').removeEmpty().getPhysio('puls');   puls = puls(:);   missing_puls = cellfun( 'isempty', puls(:).getPath() );
resp = e.getSerie('phy').removeEmpty().getPhysio('resp');   resp = resp(:);   missing_resp = cellfun( 'isempty', resp(:).getPath() );

run_all = e.getSerie('fmri_nm_bold').removeEmpty();

idx_missing = missing_info | missing_puls | missing_resp;
%%%%% !!!!
% if ~any(idx_missing) % only good complete data
%     idx_missing = logical(size(missing_info));
% end
%%%%% !!!!

idx_ok = ~idx_missing;

run_phy_missing = run_all( idx_missing );
run_phy_ok      = run_all( idx_ok );

volume_phy_ok = run_phy_ok.getVolume('^wts_OC');
outdir_phy_ok = volume_phy_ok.getDir();
rp_phy_ok     = run_phy_ok.getRP('rp_spm');
mask_phy_ok   = run_phy_ok.getExam().getSerie('anat').getVolume('^rwp[23]');
info_ok       = info( idx_ok );
puls_ok       = puls( idx_ok );
resp_ok       = resp( idx_ok );

if any(idx_missing)

    volume_phy_missing = run_phy_missing.getVolume('^wts_OC');
    outdir_phy_missing = volume_phy_missing.getDir();
    rp_phy_missing     = run_phy_missing.getRP('rp_spm');
    mask_phy_missing   = run_phy_missing.getExam().getSerie('anat').getVolume('^rwp[23]').squeeze();

end


%% Prepare job : common

cfg = global_parameters;

cfg.TR     = 1.660;
cfg.nSlice = 60;

cfg.noiseROI_thresholds   = [0.95 0.80];     % keep voxels with tissu probabilty >= 95%
cfg.noiseROI_n_voxel_crop = [2 1];           % crop n voxels in each direction, to avoid partial volume
cfg.noiseROI_n_components = 10;              % keep n PCA componenets

cfg.rp_threshold = 1.0;  % Threshold above which a stick regressor is created for corresponding volume of exceeding value

cfg.print_figures = 0; % 0 , 1 , 2 , 3

cfg.rp_order     = 24;   % can be 6, 12, 24
% 6 = just add rp, 12 = also adds first order derivatives, 24 = also adds first + second order derivatives
cfg.rp_method    = 'FD'; % 'MAXVAL' / 'FD' / 'DVARS'

cfg.display  = 0;
cfg.redo     = 0;
cfg.walltime = '04:00:00';
cfg.mem      = '4G';


%% Prepare job : ok

% ALWAYS MANDATORY
cfg.physio   = 1;
cfg.noiseROI = 1;
cfg.rp       = 1;
cfg.volume = volume_phy_ok;
cfg.outdir = outdir_phy_ok;

% Physio
cfg.physio_Info = info_ok;
cfg.physio_PULS = puls_ok;
cfg.physio_RESP = resp_ok;
cfg.physio_RETROICOR        = 1;
cfg.physio_HRV              = 1;
cfg.physio_RVT              = 1;
cfg.physio_logfiles_vendor  = 'Siemens_Tics'; % Siemens CMRR multiband sequence, only this one is coded yet
cfg.physio_logfiles_align_scan = 'last';         % 'last' / 'first'
% Determines which scan shall be aligned to which part of the logfile.
% Typically, aligning the last scan to the end of the logfile is beneficial, since start of logfile and scans might be shifted due to pre-scans;
cfg.physio_slice_to_realign    = 'middle';       % 'first' / 'middle' / 'last' / sliceNumber (integer)
% Slice to which regressors are temporally aligned. Typically the slice where your most important activation is expected.

% noiseROI
cfg.noiseROI_mask   = mask_phy_ok;
cfg.noiseROI_volume = volume_phy_ok;

% Realignment Parameters
cfg.rp_file = rp_phy_ok;


cfg.jobname  = 'spm_physio_ok';
job_physio_tapas( cfg );


%% Prepare job : missing

if any(idx_missing)

    % ALWAYS MANDATORY
    cfg.physio   = 0; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cfg.noiseROI = 1;
    cfg.rp       = 1;
    cfg.volume = volume_phy_missing;
    cfg.outdir = outdir_phy_missing;

    % Physio
    cfg.physio_Info = [];
    cfg.physio_PULS = [];
    cfg.physio_RESP = [];
    cfg.physio_RETROICOR = 0;
    cfg.physio_HRV       = 0;
    cfg.physio_RVT       = 0;

    % noiseROI
    cfg.noiseROI_mask   = mask_phy_missing;
    cfg.noiseROI_volume = volume_phy_missing;

    % Realignment Parameters
    cfg.rp_file = rp_phy_missing;

    cfg.jobname  = 'spm_physio_missing';
    job_physio_tapas( cfg );

end


%% END

save e.mat e
