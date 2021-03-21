function  SPM_blockdesign(scan_dir,file_type,TR,stim)


%% create output folder

newdir = [scan_dir 'SPM'];
if ~exist(newdir,'dir'); mkdir(newdir); end
addpath(newdir)

%% create spm structure including the file lists

matlabbatch{1}.spm.stats.fmri_spec.dir = {scan_dir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

matlabbatch{1}.spm.stats.fmri_spec.dir = {[newdir]};
units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.units = units;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
%%
nii = load_nii(file_type{1});
clear scans
for i = 1:size(nii.img,4)
    scans{i} = [scan_dir file_type{1} ',' num2str(i)];
end

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans';
%%
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.name = 'stim';


matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = stim.(units).onsets;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = stim.(units).duration;

matlabbatch{1}.spm.stats.fmri_spec.sess.cond.tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 100;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.length = 20;
matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.order = 4;
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

