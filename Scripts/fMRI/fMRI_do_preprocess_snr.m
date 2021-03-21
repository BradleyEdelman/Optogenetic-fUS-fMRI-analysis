function fMRI_do_preprocess_snr(func_fold,func_file,nii,anat_fold,param)

if isunix; slash = '/'; elseif ispc; slash = '\'; end
preproc = param.preproc;
smooth = param.smooth;
Dummy = param.Dummy;
order = param.order;
template = param.template;

dnii.img = nii.img(:,:,:,Dummy+1:end);
dnii = make_nii(dnii.img,[1 1 1]);
save_nii(dnii,[func_fold 'snr' func_file])

 % Realign
p1 = spm_select('ExtList', pwd, ['^' ['snr' func_file]], 1:size(dnii.img,4));
P = char(p1);
spm_realign(P);
spm_reslice(P);


% Registration to make sure func and anat are aligned
clear matlabbatch

anat_files = dir(anat_fold);
anat_files(~contains({anat_files.name},'anatomy.nii')) = [];
anat_files = anat_files(1);
ref = [anat_files.folder slash anat_files.name(1:end-4) '2_ds.nii'];

source = spm_select('ExtList', pwd, ['^' ['snr' func_file]], 1);
other = spm_select('ExtList', pwd, ['^' ['snr' func_file]], 2:size(dnii.img,4));

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {ref};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {source};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(other);
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'c';
spm_jobman('run',matlabbatch);

% Registration
ave_anat = [anat_fold 'ave_anatomy_ds.nii'];

clear matlabbatch

anat_files = dir(anat_fold);
anat_files(~contains({anat_files.name},'anatomy.nii')) = [];
anat_files = anat_files(1);
source = [anat_files.folder slash anat_files.name(1:end-4) '2_ds.nii'];

other = spm_select('ExtList', pwd, ['^' ['crad' func_file]], 1:size(dnii.img,4));

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {ave_anat};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {source};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(other);
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'c';
spm_jobman('run',matlabbatch);
        




%%
% Registration to make sure func and anat are aligned
% clear matlabbatch
% 
% anat_files = dir(anat_fold);
% anat_files(~contains({anat_files.name},'anatomy.nii')) = [];
% anat_files = anat_files(1);
% ref = [anat_files.folder slash 'ave_anatomy_ds.nii'];
% 
% source = spm_select('ExtList', pwd, ['^' ['snr' func_file]], 1);
% other = spm_select('ExtList', pwd, ['^' ['snr' func_file]], 2:size(dnii.img,4));
% 
% matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {ref};
% matlabbatch{1}.spm.spatial.coreg.estwrite.source = {source};
% matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(other);
% matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
% matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
% matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
% matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
% matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
% matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
% matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
% matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'c';
% spm_jobman('run',matlabbatch);
