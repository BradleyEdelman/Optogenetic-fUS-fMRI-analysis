function fMRI_Anat_Template(storage,sessions,idx)

if isunix; slash = '/'; elseif ispc; slash = '\'; end

% sample functional scan for size info
recon_fold = [storage sessions{1}(1:8) slash sessions{1}(17:end) slash];
func_fold = dir(recon_fold);
func_fold(~contains({func_fold.name},'f')) = [];
func_fold = func_fold(1);
func_file = [func_fold.folder slash func_fold.name slash func_fold.name(1:end-1) '_EPI.nii'];
f = load_nii(func_file);
sz = size(f.img(:,:,:,1));

% Specify reconstructed anatomical scans from session
recon_fold = [storage sessions{idx}(1:8) slash sessions{idx}(17:end) slash];
anat_fold = dir(recon_fold);
anat_fold(~contains({anat_fold.name},'a')) = [];
anat_fold = [recon_fold anat_fold(1).name slash];
anat_files = dir(anat_fold);
anat_files(~contains({anat_files.name},'anatomy.nii')) = [];
ref = {[anat_files(1).folder slash anat_files(1).name]};

% resave original and downsample for functional registration to group template
tmp = load_untouch_nii(ref{1});
tmpimg = flipud(permute(tmp.img,[2 1 3]));
ref2 = make_nii(tmpimg,[1 1 1]);
save_nii(ref2,[anat_files(1).folder slash anat_files(1).name(1:end-4) '2.nii']);

tmpimg = flipud(permute(imresize3(tmp.img,sz),[2 1 3]));
ref2 = make_nii(tmpimg,[1 1 1]);
save_nii(ref2,[anat_files(1).folder slash anat_files(1).name(1:end-4) '2_ds.nii']);



niitmp = load_untouch_nii(ref{1});
img(:,:,:,1) = niitmp.img;

refsession = sessions(idx);
sessions(idx) = [];

for i = 1:size(sessions,1)
    
    clear matlabbatch
    
    recon_fold = [storage sessions{i}(1:8) slash sessions{i}(17:end) slash];
    anat_fold = dir(recon_fold);
    anat_fold(~contains({anat_fold.name},'a')) = [];
    anat_fold = [recon_fold anat_fold(1).name slash];
    anat_files = dir(anat_fold);
    anat_files(~contains({anat_files.name},'anatomy.nii')) = [];
    source = {[anat_files(1).folder slash anat_files(1).name]};
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = ref;
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = source;
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'c';
    spm_jobman('run',matlabbatch);
    
    % registered invidual anatomical
    csource = [anat_files(1).folder slash 'c' anat_files(1).name]; 
    niitmp = load_untouch_nii(csource);
    img(:,:,:,i+1) = niitmp.img;
    
    % resave original and downsample for functional registration to group template
    tmp = load_untouch_nii(source{1});
    tmpimg = flipud(permute(tmp.img,[2 1 3]));
    source2 = make_nii(tmpimg,[1 1 1]);
    save_nii(source2,[anat_files(1).folder slash anat_files(1).name(1:end-4) '2.nii']);
    
    tmpimg = flipud(permute(imresize3(tmp.img,sz),[2 1 3]));
    source2 = make_nii(tmpimg,[1 1 1]);
    save_nii(source2,[anat_files(1).folder slash anat_files(1).name(1:end-4) '2_ds.nii']);
    
end

aveimg = flipud(permute(mean(img,4),[2 1 3])); % orig size
aveimg_ds = flipud(permute(imresize3(mean(img,4),sz),[2 1 3])); % func size
avenii = make_nii(aveimg,[1 1 1]);
avenii_ds = make_nii(aveimg_ds,[1 1 1]);

% save average anatomy in each session's anat directory
for i = 1:size(sessions,1)+1
    
    if i<=size(sessions,1)
        recon_fold = [storage sessions{i}(1:8) slash sessions{i}(17:end) slash];
    else
        recon_fold = [storage refsession{1}(1:8) slash refsession{1}(17:end) slash]; 
    end
        
    anat_fold = dir(recon_fold);
    anat_fold(~contains({anat_fold.name},'a')) = [];
    anat_fold = [recon_fold anat_fold(1).name slash];
    anat_files = dir(anat_fold);
    anat_files(~contains({anat_files.name},'anatomy.nii')) = [];
    ave_file = [anat_files(1).folder slash 'ave_anatomy.nii'];
    ave_file_ds = [anat_files(1).folder slash 'ave_anatomy_ds.nii'];
    
    save_nii(avenii,ave_file);
    save_nii(avenii_ds,ave_file_ds);
    
end
