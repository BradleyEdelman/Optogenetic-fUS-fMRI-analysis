function SPM_blockcontrast(scan_dir)

spm_fold = [scan_dir 'SPM/'];
cd(spm_fold);
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spm_fold 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method = struct('Classical',1);      
spm_jobman('initcfg');
spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);


%%
clear matlabbatch

matlabbatch{1}.spm.stats.con.spmmat = {[spm_fold 'SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 1; % delete existing contrasts;

load(matlabbatch{1}.spm.stats.con.spmmat{1}) % load spm file

matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_1stBasisFunc';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 0 0];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'stim_2ndBasisFunc';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [0 1 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'stim_AllBasisFunc';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [1 1 0];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';


%Run it
spm_jobman('run',matlabbatch);
