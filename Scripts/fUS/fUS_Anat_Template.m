function fUS_Anat_Template(storage,sessions,idx)

% Main data directory
if isunix
    Dir = '/media/bradley/Seagate Backup Plus Drive/';
    slash = '/';
elseif ispc
    Dir = 'D:\';
    slash = '\';
end
if isunix; slash = '/'; elseif ispc; slash = '\'; end

dataFold = [Dir 'FUSi' slash sessions{idx} slash];
anat_data = dir(dataFold);
anat_data(~contains({anat_data.name},'IQdata')) = [];
anat_data(contains({anat_data.name},'P')) = [];
anatFile = [anat_data(end).folder slash anat_data(end).name];
[adata,apdi,apsd,aiqr] = Clutterfilt_preproc(anatFile,'method','nosvd','filtcutoff',30);
ref = 10*log10(apdi./max(max(apdi)));

% all respective flow analysis formats
flow(1) = fUS_FlowAnalysis(aiqr,ref);
flow_ds(1).veloc = flow(1).veloc(1:2:end,1:2:end);
flow_ds(1).I_neg = flow(1).I_neg(1:2:end,1:2:end);
flow_ds(1).I_pos = flow(1).I_pos(1:2:end,1:2:end);
flow_ds(1).signed = flow(1).signed(1:2:end,1:2:end);

flowc(1) = flow(1);
flowc_ds(1).veloc = flowc(1).veloc(1:2:end,1:2:end);
flowc_ds(1).I_neg = flowc(1).I_neg(1:2:end,1:2:end);
flowc_ds(1).I_pos = flowc(1).I_pos(1:2:end,1:2:end);
flowc_ds(1).signed = flowc(1).signed(1:2:end,1:2:end);

% store original and registered images for all animals
img(:,:,1) = ref;
cimg(:,:,1) = ref;

refsession = sessions(idx);
sessions(idx) = [];

for i = 1:size(sessions,1)
    
    recon_fold = [storage sessions{i}(1:8) slash sessions{i}(10:end) slash];
    if ~exist(recon_fold,'dir'); mkdir(recon_fold); end
    
    dataFold = [Dir 'FUSi' slash sessions{i} slash];
    anat_data = dir(dataFold);
    anat_data(~contains({anat_data.name},'IQdata')) = [];
    anat_data(contains({anat_data.name},'P')) = [];
    anatFile = [anat_data(end).folder slash anat_data(end).name];
    [adata,apdi,apsd,aiqr] = Clutterfilt_preproc(anatFile,'method','svd','filtcutoff',20);
    source = 10*log10(apdi./max(max(apdi)));
    flow(i+1) = fUS_FlowAnalysis(aiqr,source);
    flow_ds(i+1).veloc = flow(i+1).veloc(1:2:end,1:2:end);
    flow_ds(i+1).I_neg = flow(i+1).I_neg(1:2:end,1:2:end);
    flow_ds(i+1).I_pos = flow(i+1).I_pos(1:2:end,1:2:end);
    flow_ds(i+1).signed = flow(i+1).signed(1:2:end,1:2:end);
    
    % Format flow images for applying registration transformation
    other(:,:,1) = flow(i+1).veloc;
    other(:,:,2) = flow(i+1).I_neg;
    other(:,:,3) = flow(i+1).I_pos;
    other(:,:,4) = flow(i+1).signed;
    
    % registered invidual anatomical
    [ref,source,csource,other,otherc] = fUS_Coreg(ref,source,other);
    save([recon_fold 'cAnat.mat'],'csource','-v7.3')
    img(:,:,i+1) = source;
    cimg(:,:,i+1) = csource;
    flowc(i+1).veloc = otherc(:,:,1); flowc_ds(i+1).veloc = flowc(i+1).veloc(1:2:end,1:2:end);
    flowc(i+1).I_neg = otherc(:,:,2); flowc_ds(i+1).I_neg = flowc(i+1).I_neg(1:2:end,1:2:end);
    flowc(i+1).I_pos = otherc(:,:,3); flowc_ds(i+1).I_pos = flowc(i+1).I_pos(1:2:end,1:2:end);
    flowc(i+1).signed = otherc(:,:,4); flowc_ds(i+1).signed = flowc(i+1).signed(1:2:end,1:2:end);
    
    % resave original and downsample for functional registration to group template
    save([recon_fold 'Anat.mat'],'source','-v7.3')
    source2 = source(1:2:end,1:2:end);
    save([recon_fold 'Anat_ds.mat'],'source2','-v7.3')
    
    % save flow analysis for original and registered
    F = flow(i+1); Fc = flowc(i+1); Fds = flow_ds(i+1); Fcds = flowc_ds(i+1);
    save([recon_fold 'Flow.mat'],'F','-v7.3');
    save([recon_fold 'Flowc.mat'],'Fc','-v7.3');
    save([recon_fold 'Flow_ds.mat'],'Fds','-v7.3');
    save([recon_fold 'Flowc_ds.mat'],'Fcds','-v7.3');
    
end

% also save images to reference anatomy folder
recon_fold = [storage refsession{1}(1:8) slash refsession{1}(10:end) slash];
if ~exist(recon_fold,'dir'); mkdir(recon_fold); end
save([recon_fold 'cAnat.mat'],'ref','-v7.3')
save([recon_fold 'Anat.mat'],'ref','-v7.3')
source2 = ref(1:2:end,1:2:end);
save([recon_fold 'Anat_ds.mat'],'source2','-v7.3')
% and flow analysis
F = flow(1); Fc = flowc(1); Fds = flow_ds(1); Fcds = flowc_ds(1);
save([recon_fold 'Flow.mat'],'F','-v7.3');
save([recon_fold 'Flowc.mat'],'Fc','-v7.3');
save([recon_fold 'Flow_ds.mat'],'Fds','-v7.3');
save([recon_fold 'Flowc_ds.mat'],'Fcds','-v7.3');


aveimg = mean(cimg,3);
aveimg_ds = aveimg(1:2:end,1:2:end);
totimg = cimg;
totimg_ds = cimg(1:2:end,1:2:end,:);
totimg_orig = img;
totimg_orig_ds = img(1:2:end,1:2:end,:);

for i = 1:size(sessions,1)+1
    
    if i <= size(sessions,1)
    	recon_fold = [storage sessions{i}(1:8) slash sessions{i}(10:end) slash];
    else
        recon_fold = [storage refsession{1}(1:8) slash refsession{1}(10:end) slash];
    end
    
    save([recon_fold 'Ave_Anat.mat'],'aveimg','-v7.3')
    save([recon_fold 'Ave_Anat_ds.mat'],'aveimg_ds','-v7.3')
    
    % also save all images to register to a particular brain
    save([recon_fold 'Total_Anat.mat'],'totimg','-v7.3')
    save([recon_fold 'Total_Anat_ds.mat'],'totimg_ds','-v7.3')
    
    % and original images
    save([recon_fold 'Total_Anat_orig.mat'],'totimg_orig','-v7.3')
    save([recon_fold 'Total_Anat_orig_ds.mat'],'totimg_orig_ds','-v7.3')
    
end
