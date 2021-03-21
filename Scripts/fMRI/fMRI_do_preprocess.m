function FUNC = fMRI_do_preprocess(func_fold,func_file,nii,anat_fold,param)

if isunix; slash = '/'; elseif ispc; slash = '\'; end
preproc = param.preproc;
smooth = param.smooth;
Dummy = param.Dummy;
order = param.order;
template = param.template;

if isequal(preproc,1)
    
    dnii.img = nii.img(:,:,:,Dummy+1:end);
    dnii = make_nii(dnii.img,[1 1 1]);
    save_nii(dnii,[func_fold 'd' func_file])
    p1 = spm_select('ExtList', pwd, ['^' ['d' func_file]], 1:size(dnii.img,4));
    P = char(p1);
    
    % Slice-timing correction
    clear matlabbatch
    matlabbatch{1}.spm.temporal.st.scans = {cellstr(P)};
    matlabbatch{1}.spm.temporal.st.nslices = 24;
    matlabbatch{1}.spm.temporal.st.tr = 1.5;
    matlabbatch{1}.spm.temporal.st.ta = 1.5 - (1.5/24);
    matlabbatch{1}.spm.temporal.st.so = [0:2:22 1:2:23];
    matlabbatch{1}.spm.temporal.st.refslice = 12;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    spm_jobman('run',matlabbatch);
    
    % Realign
    p1 = spm_select('ExtList', pwd, ['^' ['ad' func_file]], 1:size(dnii.img,4));
    P = char(p1);
    spm_realign(P);
    spm_reslice(P);
    
    % Registration to make sure func and anat are aligned
    clear matlabbatch
        
    anat_files = dir(anat_fold);
    anat_files(~contains({anat_files.name},'anatomy.nii')) = [];
    anat_files = anat_files(1);
    ref = [anat_files.folder slash anat_files.name(1:end-4) '2_ds.nii'];

    source = spm_select('ExtList', pwd, ['^' ['rad' func_file]], 1);
    other = spm_select('ExtList', pwd, ['^' ['rad' func_file]], 2:size(dnii.img,4));

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
    if ~isempty(anat_fold)
        ave_anat = [anat_fold 'ave_anatomy_ds.nii'];
    else
        ave_anat = [];
    end
    
    if exist(ave_anat,'file') && isequal(template,1)
        
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
        
        % Smooth
        rfile = ['ccrad' func_file];
        srfile = SPM_smooth(rfile,smooth);

    else
        
        rfile = ['rad' func_file];
        srfile = SPM_smooth(rfile,smooth);
        
    end
    
    % Average realigned EPI image
    rnii = load_untouch_nii(rfile);
    rnii.img = uint16(mean(rnii.img,4));
    save_untouch_nii(rnii,['mean' rfile]);
    
else
    
    % Registration
    ave_anat = [anat_fold 'ave_anatomy_ds.nii'];
    if exist(ave_anat,'file')
        
        rfile = ['ccrad' func_file];
        srfile{1} = ['s' num2str(smooth) 'ccrad' func_file];
        
    else
        
        rfile = ['crad' func_file];
        srfile{1} = ['s' num2str(smooth) 'crad' func_file];
        
    end
    
end

% GLM analysis
% GLM = GLM_fMRI(func_fold,anat_fold,srfile{1},param);
GLM = GLM_fMRI(func_fold,anat_fold,rfile,param); % unsmoothed data

Data.orig = nii;
Data.realign = load_untouch_nii(rfile);
Data.smooth = load_untouch_nii(srfile{1});
Data.file = func_file;
Data.path = func_fold;
Data.saveFold = func_fold;

FUNC.type = 'fmri';
FUNC.Data = Data;
FUNC.GLM = GLM;

if isequal(template,1)
    t = '_template';
    func_fold = [func_fold 'template' slash];
else
    t = [];
    func_fold = [func_fold 'notemplate' slash];
end
if ~exist(func_fold,'dir'); mkdir(func_fold); end
funcfile = [func_fold 'Func.mat'];
save(funcfile,'FUNC','-v7.3')



