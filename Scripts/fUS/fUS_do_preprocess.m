function FUNC = fUS_do_preprocess(func_fold,func_data,anat_fold,param)

if isunix; slash = '/'; elseif ispc; slash = '\'; end
preproc = param.preproc;
smooth = param.smooth;
Dummy = param.Dummy;
order = param.order;
template = param.template;
templateidx = param.templateidx;
shift = param.shift;

if isequal(preproc,1)
    
    % Clutter Filter
    if isstruct(func_data)
        funcFile = [func_data.folder slash func_data.name];
        [fdata,fpdi,fpsd] = Clutterfilt_preproc(funcFile,'method','nosvd','filtcutoff',30);
    else
        fdata = [];
        fpdi = func_data;
        fpsd = [];
    end
    
    % Adjust for FOV shift error (hard coded in 'Preprocessing_fUS.m')
    shift = param.shift;
    for i = 1:size(fpdi,3)
        fpdi(:,:,i) = circshift(fpdi(:,:,i),shift);
    end
    
    % Registration
    if ~isempty(anat_fold)
        if isnan(templateidx)
            ave_anat = [anat_fold 'Ave_Anat_ds.mat'];
        else
            ave_anat = [anat_fold 'Total_Anat_orig_ds.mat'];
        end
    else
        ave_anat = [];
    end
    
    if exist(ave_anat,'file') && isequal(template,1)
        
        % load anat
        tmp = load(ave_anat);
        if isnan(templateidx)
            ref = tmp.aveimg_ds;
        else
            ref = tmp.totimg_orig_ds(:,:,templateidx);
        end
        
        % load func
        anat_files = dir(anat_fold);
        anat_files(~contains({anat_files.name},'Anat')) = [];
        anat_files = anat_files(1);
        source = [anat_files.folder slash anat_files.name(1:end-4) '_ds.mat'];
        tmp = load(source); source = tmp.source2;
        
        % ensure anat and func are same size
        fpdi2 = imresize3(fpdi,[size(ref) size(fpdi,3)]);
        other = fpdi2;
        
        % register animal anat and func to selected reference
        [ref,source,csource,fpdi3,cfpdi3] = fUS_Coreg(ref,source,other);
        
        % Average low res PDI image
        ave_cfpdi3 = mean(cfpdi3,3);
        if isstruct(func_data)
            save([func_fold func_data.name(end-16:end-9) '_Ave_Func.mat'],'ave_cfpdi3','-v7.3')
        else
            data = dir(func_fold);
            data(~contains({data.name},'Func')) = [];
            data = data(end);
            save([func_fold data.name(1:3) '_Ave_Func.mat'],'ave_cfpdi3','-v7.3')
        end
    else
        
        cfpdi3 = fpdi;
        
    end
    
    % smooth
    for i = 1:240
        cfpdi3(:,:,i) = medfilt2(cfpdi3(:,:,i),[3 3]);
    end
    
    
    % GLM analysis
    [Block,GLM] = GLM_fUS(func_fold,anat_fold,cfpdi3,param);    
        
        
        
end

Data.pdi = cfpdi3;
Data.pdi_orig = fpdi;
Data.psd = fpsd;
if isstruct(func_data)
    Data.file = func_data.name;
    Data.path = func_data.folder;
else
    data = dir(func_fold);
    data(~contains({data.name},'Func')) = [];
    data = data(end);
    
    Data.file = data.name;
    Data.path = data.folder;
end

FUNC.type = 'fus';
FUNC.Data = Data;
FUNC.Block = Block;
FUNC.GLM = GLM;

if isequal(template,1)
    func_fold = [func_fold 'template' slash];
else
    func_fold = [func_fold 'notemplate' slash];
end
if ~exist(func_fold,'dir'); mkdir(func_fold); end

if isstruct(func_data)
    func_file_save = [func_fold func_data.name(end-16:end-9) '_Func.mat'];
    func_file_iqr = [func_fold func_data.name(end-16:end-9) '_IQR.mat'];
else
    func_file_save = [func_fold data.name(1:3) '_Func.mat'];
    func_file_iqr = [func_fold data.name(1:3) '_IQR.mat'];
end
    
save(func_file_save,'FUNC','-v7.3')

% IQR = fdata.D.IQR;
% save(iqr_file_save,'IQR','-v7.3')


