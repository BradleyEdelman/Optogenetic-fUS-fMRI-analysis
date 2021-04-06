%% Bulk preprocessing script

% Main data directory
if isunix
    Dir = '/media/bradley/Seagate Backup Plus Drive/';
    slash = '/';
elseif ispc
    Dir = 'D:\';
    slash = '\';
end
% addpath(genpath([Dir 'Preprocessing' slash]));
storage = [Dir 'Data_Processed' slash 'fUS' slash ]; if ~exist(storage,'dir'); mkdir(storage); end
addpath(genpath(storage));

% Thy1
base_fold = {'20191119_4346075_N_D2';
    '20191120_4364143_R1_D2';
    '20191120_4364124_L1_D2';
    '20191123_4364123_R1_D3';
    '20191123_4364122_N_D3';
    '20200117_4389772_N_D2';
    '20200117_4389773_R1_D2';
    };

shift = [-3 0; 0 0; -1 0; -2 0; 0 0; 0 0; 0 0];

SNR = 1;
x1 = {'0.1'; {'20-11-48' '20-24-46'}; {'10-39-27' '10-51-23'}; {'19-06-27' '19-21-23'}; {'11-31-43' '11-47-42'}; {'14-28-40' '14-40-51'}; {'10-10-25' '10-24-04'}; {'17-20-36'}};
x2 = {'0.5'; {'19-11-52' '19-32-40'}; {'11-02-58' '11-15-05'}; {'20-00-17' '20-15-10'}; {'12-23-14'}; {'15-05-35' '15-19-32'}; {'10-36-32'}; {'17-08-15'}}; 
x3 = {'1.0'; {'19-46-05' '19-58-55'}; {'11-27-12' '11-38-48'}; {'19-34-59' '19-47-22'}; {'11-59-30' '12-11-10'}; {'14-01-17' '14-15-52'}; {'10-49-11'}; {'16-55-19'}};
stimidx = {x1, x2, x3};

% % Control
% base_fold = {'20200224_4419409_N';
%     '20200224_4419410_N';
%     };
% 
% shift = [0 0; 0 0];
% 
% SNR = 1;
% x1 = {'0.1'; {'14-39-08' '14-51-23'}; {'17-27-07' '18-03-47'}};
% x2 = {'0.5'; {'14-15-09' '14-27-30'}; {'18-18-38' '18-31-45'}};
% x3 = {'1.0'; {'15-03-03' '15-14-46'}; {'17-27-07' '17-39-26'}};
% stimidx = {x1, x2, x3};


% 
idx = 1;
% fUS_Anat_Template(storage,base_fold,idx)

%%

for i = 1:size(base_fold,1)
    
    dataFold = [Dir 'FUSi' slash base_fold{i} slash];
    
    func_data = dir(dataFold);
    func_data(~contains({func_data.name},'Func')) = [];
    anat_data = dir(dataFold);
    anat_data(~contains({anat_data.name},'IQdata')) = [];
    anat_data(contains({anat_data.name},'P')) = [];
    
    cd(dataFold)
    
    recon_fold = [storage base_fold{i}(1:8) slash base_fold{i}(10:end) slash];
    if ~exist(recon_fold,'dir'); mkdir(recon_fold); end
        
    param.preproc = 1;
    param.smooth = 2.5; % um *100
    param.Dummy = 10; %# TRs
    param.order = 4; %GLM gamma
    param.template = 1; % Use study brain template
    param.templateidx = 1; % Use a single brain (number), or average (nan)
    param.shift = shift(i,:); % FOV correction
    
%     anatFile = [anat_data(end).folder slash anat_data(end).name];
%     [adata,apdi,apsd] = Clutterfilt_preproc(anatFile,'method','nosvd','filtcutoff',30);
%     apdi = 10*log10(apdi./max(max(apdi)));
%     
    if isequal(SNR,1) % Average runs for higher SNR
            
        func_opt = {func_data.name};
        for k = 1:size(stimidx,2) % stim intensities
            
            snr_fold = [recon_fold stimidx{k}{1} slash];
            if ~exist(snr_fold,'dir'); mkdir(snr_fold); end
            func_fold = snr_fold;
            func_file = [snr_fold stimidx{k}{1} '_Func.mat'];
            
            clear fdata fpdi fpsd
            if ~exist(func_file,'file')
            
                snrfiles = stimidx{k}{i+1};
                snridx = contains(func_opt,snrfiles);
                snridx = find(snridx == 1);
                for l = 1%:size(snridx,2)
                    snr_file = [dataFold func_opt{snridx(l)}];
                    [fdata{l},fpdi(:,:,:,l),fpsd(:,:,l)] = Clutterfilt_preproc(snr_file,'method','nosvd','filtcutoff',30);
                end

                cd(snr_fold)
                func_data = mean(fpdi,4);

                save(func_file,'fpdi','-v7.3');
            else
                
                load(func_file)
                func_data = mean(fpdi,4);
                
            end
            
            anat_fold = recon_fold;
            
            FUNC = fUS_do_preprocess(func_fold,func_data,anat_fold,param);
            
        end
        
    else
        
        for j = size(func_data,1)
            
            func_fold = [recon_fold func_data(j).name(end-16:end-9) slash];
            if ~exist(func_fold,'dir'); mkdir(func_fold); end
            cd(func_fold)
            
            anat_fold = recon_fold;
            
%             func_file = [func_fold 'template\' func_data(j).name(end-16:end-9) '_Func.mat'];
%             if exist(func_file,'file')
%                 load(func_file)
%             else
                FUNC = fUS_do_preprocess(func_fold,func_data(j),anat_fold,param);
%             end
            
%             [S{i},T{i},R{i}] = fusi_compute_snr(FUNC.Data.pdi);
                    
        end
    end
%     pause
end

% save('D:\Data_Processed\fUSI_SNR.mat','S','T','R')
%% Individual and group level fixed effects, time series

clear all
% close all

% Main data directory
if isunix
    Dir = '/media/bradley/Seagate Backup Plus Drive/';
    slash = '/';
elseif ispc
    Dir = 'D:\';
    slash = '\';
end
% addpath(genpath([Dir 'Preprocessing' slash]));

addpath(genpath('C:\Users\bedelman\Documents\GitHub\ofUSI\'))
storage = [Dir 'Data_Processed' slash 'fUS' slash ]; if ~exist(storage,'dir'); mkdir(storage); end
addpath(genpath(storage));

% Thy1
base_fold = {'20191119_4346075_N_D2';
    '20191123_4364122_N_D3';
    '20191123_4364123_R1_D3';
    '20191120_4364124_L1_D2';
    '20191120_4364143_R1_D2';
    '20200117_4389772_N_D2';
    '20200117_4389773_R1_D2';
    };
descr = 'Thy1';

% Individual
% '20191120_4364124_L1_D2' Vascular M1
% base_fold = {'20191119_4346075_N_D2'; %'20191120_4364124_L1_D2','20200117_4389772_N_D2'
%     '20191120_4364124_L1_D2';
%     '20200117_4389772_N_D2';
%     '20200117_4389773_R1_D2';
%     };

% '20191123_4364122_N_D3'; no
% '20191120_4364124_L1_D2'; yes
% '20191120_4364143_R1_D2'; no
% '20191119_4346075_N_D2'; yes
% '20191123_4364123_R1_D3'; no
% base_fold = {'20191123_4364123_R1_D3';}
% descr = 'single';

% Control
% base_fold = {'20200224_4419409_N';
%     '20200224_4419410_N';
%     };
% descr = 'ctrl';

% Ephys Subset
% base_fold = {'20191119_4346075_N_D2';
%     '20191123_4364123_R1_D3';
%     '20191123_4364122_N_D3';
%     '20191120_4364143_R1_D2';
%     };
% descr = 'ephys';

%%
fusroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2'...
    'RCg1' 'RCg2' 'RM2' 'LS1BF' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

param.order = 4; %GLM gamma
param.template = 1; % User study brain template; must have template for group level
param.Dummy = 10; %# TRs
param.Pthresh = 0.005;
param.maxt = 40;
param.fusroi = fusroi;
param.descr = descr;
SNR = 1; % Always stim intensity-based analysis
%%
% Load and compile data from individual animals
fus_organize_data(storage,base_fold,slash,param)
% Fixed effects analysis for each animal
fus_individual_fixed_effects(storage,base_fold,slash,param)
% Sig active voxel count (+ flow count)
fus_active_voxel_count(storage,base_fold,slash,param)
% Fixed effects analysis at the group level
fus_group_fixed_effects(storage,base_fold,slash,param)
% Plot group activation maps
fus_plot_group_fixed_effects(storage,base_fold,slash,param)
% Extract time series for roi and vascular components thereof
fus_timeseries_analysis(storage,base_fold,slash,param)
% Plot time series info
fus_plot_timeseries_analysis(storage,base_fold,slash,param)
%%
% Plot ROI activation (percent, t-stat)
fus_plot_ROI_info(storage,base_fold,slash,param)
% Regression Analysis on time series info
fus_regression_analysis(storage,base_fold,slash,param)
% Regression Analysis on time series info for flow components
fus_regression_analysis_flow(storage,base_fold,slash,param)
%%

fusroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2'...
    'RCg1' 'RCg2' 'RM2' 'LS1BF' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};
fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
    'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};


tfmri = load('D:\Data_Processed\fMRI\pre_regression.mat'); tfmri = tfmri.B_roi_sig;
tfus = load('D:\Data_Processed\fUS\Thy1_regression.mat'); tfus = tfus.B_roi_sig;
figure; hold on
for i = 1:20 
    
    ref = fusroi{i};
    fmriidx = find(strcmp(fmriroi,ref));
    
    scatter(tfmri(fmriidx),tfus(i));
    tfmritmp(i) = tfmri(fmriidx); tfustmp(i) = tfus(i);
    text(tfmri(fmriidx)+0.1,tfus(i),fusroi{i})
    plot(-1:10,-1:10)
    
end
xlabel('fMRI'); ylabel('fUSI');
title('ROI Sensitivity: fUSI vs fMRI');

empt = [];
empt = find(tfustmp == 0); empt = [empt find(tfmritmp == 0)];
empt = unique(empt);

Rat = abs(tfustmp)./abs(tfmritmp); Rat([6,20,empt]) = [];
fusroitmp = fusroi; fusroitmp([6,20,empt]) = [];
[B,I] = sort(Rat);
figure; barh(Rat(I)-1); grid on
set(gca,'ytick',1:size(fusroitmp,2),'yticklabel',fusroitmp(I),'xlim',[-1 2]);
title('Beta Ratio')


tfmri = load('D:\Data_Processed\fMRI\pre_regression.mat'); tfmri = tfmri.B_roi_int_sig;
tfus = load('D:\Data_Processed\fUS\Thy1_regression.mat'); tfus = tfus.B_roi_int_sig;
figure; hold on
for i = 1:20 
    
    ref = fusroi{i};
    fmriidx = find(strcmp(fmriroi,ref));
    
    scatter(tfmri(fmriidx),tfus(i));
    tfmritmp(i) = tfmri(fmriidx); tfustmp(i) = tfus(i);
    text(tfmri(fmriidx)+0.1,tfus(i),fusroi{i})
    plot(-1:3.5,-1:3.5)
    
end
xlabel('fMRI'); ylabel('fUSI');
title('ROI Intercept: fUSI vs fMRI');

empt = [];
empt = find(tfustmp == 0); empt = [empt find(tfmritmp == 0)];
empt = unique(empt);

Rat = tfustmp - tfmritmp; Rat([6,20,empt]) = [];
fusroitmp = fusroi; fusroitmp([6,20,empt]) = [];
[B,I] = sort(Rat);
figure; barh(Rat(I)); grid on
set(gca,'ytick',1:size(fusroitmp,2),'yticklabel',fusroitmp(I));
title('Intercept Ratio')

% Pixel beta and intercept
tfmri = load('D:\Data_Processed\fMRI\pre_regression.mat');
tfus = load('D:\Data_Processed\fUS\Thy1_regression.mat');

figure;
x = [tfmri.B_pix_sig(:); tfus.B_pix_sig(:)];
g = [repmat({'fmri'},size(tfmri.B_pix_sig(:),1),1); repmat({'fusi'},size(tfus.B_pix_sig(:),1),1)];
subplot(1,2,1); boxplot(x,g); title('Beta')
[h,p] = ranksum(tfmri.B_pix_sig(:),tfus.B_pix_sig(:))

x = [tfmri.B_int_pix_sig(:); tfus.B_int_pix_sig(:)];
g = [repmat({'fmri'},size(tfmri.B_int_pix_sig(:),1),1); repmat({'fusi'},size(tfus.B_int_pix_sig(:),1),1)];
subplot(1,2,2); boxplot(x,g)  ; title('Int') 
[h,p] = ranksum(tfmri.B_int_pix_sig(:),tfus.B_int_pix_sig(:))

%%
% SNR
fus = load('D:\Data_Processed\fUSI_SNR.mat');
fmri = load('D:\Data_Processed\fMRI_SNR.mat');

fustmp = vertcat(fus.T{:});
fmritmp = vertcat(fmri.T{:});

for i = 1:20
    
    ref = fus.R{1}{i};
    fmriidx = find(strcmp(fmri.R{1},ref));
    
    Tmean(i,1) = mean(fustmp(:,i));
    Tstd(i,1) = std(fustmp(:,i));
    
    Tmean(i,2) = mean(fmritmp(:,fmriidx));
    Tstd(i,2) = std(fmritmp(:,fmriidx));
    
end

%%
% 2nd-level analysis - effect of stimulation

% REGRESSION
condition = repmat([0.1 0.5 1.0]',size(base_fold,1),1);
condition = condition - repmat(mean(condition),size(condition,1),1);
reg_D = zeros(3*size(base_fold,1),size(base_fold,1));
for i = 1:size(base_fold,1)
    reg_D(1 + 3*(i-1):3 + 3*(i-1),i) = 1;
end
reg_D = [condition reg_D];

% Set up data using individual animal contrast images
reg_Data = [];
for i = 1:size(base_fold,1)
    for j = 1:size(stim,2)
        reg_Data = [reg_Data; reshape(ind_cont(:,:,i,j),1,size(ind_cont,1)*size(ind_cont,2))];
    end
end


reg_c = [1 zeros(1,size(base_fold,1))]'; % Regression contrast
reg_t_stat = zeros(1,size(reg_Data,2));
reg_p2tail = zeros(1,size(reg_Data,2));
reg_beta = zeros(size(reg_Data,2),size(reg_c,1));
for i = 1:size(reg_Data,2)
    Y = reg_Data(:,i)';

    DF = size(reg_Data,1) - 2;

    beta_hat=inv(reg_D'*reg_D)*reg_D'*Y';
    reg_beta(i,:) = reshape(beta_hat,1,size(reg_c,1));
    Var_e=(Y'-reg_D*beta_hat)'*(Y'-reg_D*beta_hat)/DF;

    %Hypothesis testing; Compute the t statistic
    reg_t_stat(i)=reg_c'*beta_hat/sqrt(Var_e*reg_c'*inv(reg_D'*reg_D)*reg_c);
    
    if isnan(reg_t_stat(i))
        reg_p2tail(i) = 1;
    else
        reg_p2tail(i) = 1 - tdist2T(reg_t_stat(i),DF);
    end
    
end
reg_t_stat = reshape(reg_t_stat,size(ind_cont,1),size(ind_cont,2));
reg_p2tail = reshape(reg_p2tail,size(ind_cont,1),size(ind_cont,2));

reg_c_p2tail = reg_p2tail;

maxt = 10;
plotmin = -100; plotmax = 0;
Cbar = [gray(abs(plotmax - plotmin)*100); fireice(2*maxt*100)]; CAX = [plotmin - maxt maxt];
Cbarp = [jet(round(1.1*100));[0 0 0]];

atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fUS/20191119/4346075_N_D2/Anat.mat';
anii = load(atmp); anii = anii.ref;
mask = ones(size(anii));
% Resize time-series param to anat sampling
reg_p2tail_tmp = imresize(reg_c_p2tail,[size(anii,1) size(anii,2)]); reg_p2tail_tmp(mask == 0) = 1.1;
reg_t_stat_tmp = imresize(reg_t_stat,[size(anii,1) size(anii,2)]); ttmp(mask == 0) = 0;

% rescale dynamic range of anat for vis
Pthresh = 0.005;
anat = anii; anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
anat = anat - maxt;
reg_anat = anat; reg_anat(reg_p2tail_tmp < Pthresh) = reg_t_stat_tmp(reg_p2tail_tmp < Pthresh); reg_anat(mask == 0) = 0;

figure;
subplot(2,2,1); imagesc(reg_D); set(gca,'yticklabel',[],'xticklabel',[]); title('Regression Design Matrix'); colormap(gca,'jet')
subplot(2,2,2); imagesc(reg_t_stat_tmp); set(gca,'yticklabel',[],'xticklabel',[]); title('t map'); colormap(gca,jet); caxis([-maxt maxt])
subplot(2,2,3); imagesc(reg_p2tail_tmp); set(gca,'yticklabel',[],'xticklabel',[]); title('p-val 2 tail'); colormap(gca,Cbarp); caxis([0 1])
subplot(2,2,4); imagesc(reg_anat); set(gca,'yticklabel',[],'xticklabel',[]); title('t map 2 tail'); colormap(gca,Cbar); caxis(CAX)


% ANOVA + POST HOC T-TESTS
condition = zeros(size(base_fold,1)*3,3);
condition(1:3:end,1) = 1; condition(2:3:end,2) = 1; condition(3:3:end,3) = 1;
anov_D = [];
for i = 1:size(base_fold,1)
    anov_D = blkdiag(anov_D,ones(3,1));
end
anov_D = [condition anov_D];

% Set up data using individual animal contrast images
anov_DATA = [];
for i = 1:size(base_fold,1)
    for j = 1:size(stim,2)
        anov_DATA = [anov_DATA; reshape(ind_cont(:,:,i,j),1,size(ind_cont,1)*size(ind_cont,2))];
    end
end

anov_c = [1 0 -1 zeros(1,size(base_fold,1)); 0 1 -1 zeros(1,size(base_fold,1))]'; % ANOVA contrast
ttest1_c = [-1 1 0 zeros(1,size(base_fold,1))]'; % t-test contrasts
ttest2_c = [0 -1 1 zeros(1,size(base_fold,1))]';
ttest3_c = [-1 0 1 zeros(1,size(base_fold,1))]';
anov_beta = zeros(size(anov_DATA,2),size(anov_c,1));
anov_F = zeros(1,size(anov_DATA,2));
anov_p = zeros(1,size(anov_DATA,2));
ttest1_t_stat = zeros(1,size(anov_DATA,2)); ttest1_p = zeros(1,size(anov_DATA,2));
ttest2_t_stat = zeros(1,size(anov_DATA,2)); ttest2_p = zeros(1,size(anov_DATA,2));
ttest3_t_stat = zeros(1,size(anov_DATA,2)); ttest3_p = zeros(1,size(anov_DATA,2));
for i = 1:size(anov_DATA,2)
    
    Y = anov_DATA(:,i)';

    DF = size(anov_D,1) - 2;

    beta_hat=pinv(anov_D'*anov_D)*anov_D'*Y';
    r = rank(anov_c);
    anov_beta(i,:) = reshape(beta_hat,1,size(anov_c,1));
    Var_e=(Y'-anov_D*beta_hat)'*(Y'-anov_D*beta_hat)/DF;
    
    anov_F(i) = (anov_c'*beta_hat)'/(r*Var_e*anov_c'*pinv(anov_D'*anov_D)*anov_c)*(anov_c'*beta_hat);
    anov_p(i) = 1 - fcdf(anov_F(i),DF,r);
    
    ttest1_t_stat(i) = ttest1_c'*beta_hat/sqrt(Var_e*ttest1_c'*pinv(anov_D'*anov_D)*ttest1_c);
    if isnan(ttest1_t_stat(i)); ttest1_p(i) = 1; else; ttest1_p(i) = 1 - tdist2T(ttest1_t_stat(i),DF); end
    
    ttest2_t_stat(i) = ttest2_c'*beta_hat/sqrt(Var_e*ttest2_c'*pinv(anov_D'*anov_D)*ttest2_c);
    if isnan(ttest2_t_stat(i)); ttest2_p(i) = 1; else; ttest2_p(i) = 1 - tdist2T(ttest2_t_stat(i),DF); end
    
    ttest3_t_stat(i) = ttest3_c'*beta_hat/sqrt(Var_e*ttest3_c'*pinv(anov_D'*anov_D)*ttest3_c);
    if isnan(ttest3_t_stat(i)); ttest3_p(i) = 1; else; ttest3_p(i) = 1 - tdist2T(ttest3_t_stat(i),DF); end
    
end
anov_F = reshape(anov_F,size(ind_cont,1),size(ind_cont,2));
anov_p = reshape(anov_p,size(ind_cont,1),size(ind_cont,2));
ttest1_t_stat = reshape(ttest1_t_stat,size(ind_cont,1),size(ind_cont,2));
ttest2_t_stat = reshape(ttest2_t_stat,size(ind_cont,1),size(ind_cont,2));
ttest3_t_stat = reshape(ttest3_t_stat,size(ind_cont,1),size(ind_cont,2));
ttest1_p = reshape(ttest1_p,size(ind_cont,1),size(ind_cont,2));
ttest2_p = reshape(ttest2_p,size(ind_cont,1),size(ind_cont,2));
ttest3_p = reshape(ttest3_p,size(ind_cont,1),size(ind_cont,2));

maxt = 10; maxf = 45;
plotmin = -100; plotmax = 0;
Cbart = [gray(abs(plotmax - plotmin)*100); fireice(2*maxt*100)]; CAXt = [plotmin - maxt maxt];
Cbarf = [gray(abs(plotmax - plotmin)*100); fireice(2*maxf*100)]; CAXf = [plotmin - maxf maxf];
Cbarp = [jet(round(1.1*100));[0 0 0]];

atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fUS/20191119/4346075_N_D2/Anat.mat';
anii = load(atmp); anii = anii.ref;
mask = ones(size(anii));
% Resize time-series param to anat sampling
anov_F_tmp = imresize(anov_F,[size(anii,1) size(anii,2)]);
anov_p_tmp = imresize(anov_p,[size(anii,1) size(anii,2)]);
ttest1_p_tmp = imresize(ttest1_p,[size(anii,1) size(anii,2)]);
ttest1_t_stat_tmp = imresize(ttest1_t_stat,[size(anii,1) size(anii,2)]);
ttest2_p_tmp = imresize(ttest2_p,[size(anii,1) size(anii,2)]);
ttest2_t_stat_tmp = imresize(ttest2_t_stat,[size(anii,1) size(anii,2)]);
ttest3_p_tmp = imresize(ttest3_p,[size(anii,1) size(anii,2)]);
ttest3_t_stat_tmp = imresize(ttest3_t_stat,[size(anii,1) size(anii,2)]);

% rescale dynamic range of anat for vis
Pthresh = 0.005;
anat = anii; anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
F_anat = anat - maxf;
F_anat(anov_p_tmp < Pthresh) = anov_F_tmp(anov_p_tmp < Pthresh);
anatt = anat - maxt;
ttest1_anat = anatt; ttest1_anat(ttest1_p_tmp < Pthresh) = ttest1_t_stat_tmp(ttest1_p_tmp < Pthresh);
ttest2_anat = anatt; ttest2_anat(ttest2_p_tmp < Pthresh) = ttest2_t_stat_tmp(ttest2_p_tmp < Pthresh);
ttest3_anat = anatt; ttest3_anat(ttest3_p_tmp < Pthresh) = ttest3_t_stat_tmp(ttest3_p_tmp < Pthresh);

figure;
subplot(5,2,1); imagesc(anov_D); set(gca,'xticklabel',[],'yticklabel',[]); title('ANOVA Design Matrix'); colormap(gca,jet);
subplot(5,2,3); imagesc(anov_p_tmp); set(gca,'xticklabel',[],'yticklabel',[]); title('F-test (ANOVA) p-val'); colormap(gca,Cbarp); caxis([0 1.1])
subplot(5,2,4); imagesc(F_anat); set(gca,'xticklabel',[],'yticklabel',[]); title('F-val (ANOVA)'); colormap(gca,Cbarf); caxis(CAXf)
subplot(5,2,5); imagesc(ttest1_p_tmp); set(gca,'xticklabel',[],'yticklabel',[]); title('p-val 0.5 - 0.1 mW'); colormap(gca,Cbarp); caxis([0 1.1])
subplot(5,2,6); imagesc(ttest1_anat); set(gca,'xticklabel',[],'yticklabel',[]); title('t-val 0.5 - 0.1 mW'); colormap(gca,Cbart); caxis(CAXt)
subplot(5,2,7); imagesc(ttest2_p_tmp); set(gca,'xticklabel',[],'yticklabel',[]); title('p-val 1.0 - 0.5 mW'); colormap(gca,Cbarp); caxis([0 1.1])
subplot(5,2,8); imagesc(ttest2_anat); set(gca,'xticklabel',[],'yticklabel',[]); title('t-val 1.0 - 0.5 mW'); colormap(gca,Cbart); caxis(CAXt)
subplot(5,2,9); imagesc(ttest3_p_tmp); set(gca,'xticklabel',[],'yticklabel',[]); title('p-val 1.0 - 0.1 mW'); colormap(gca,Cbarp); caxis([0 1.1])
subplot(5,2,10); imagesc(ttest3_anat); set(gca,'xticklabel',[],'yticklabel',[]); title('t-val 1.0 - 0.1 mW'); colormap(gca,Cbart); caxis(CAXt)



%% Time series amplitude across modality

% fusroi = {'LS2' 'LS1BF' 'LCPU' 'LS1DZ' 'LS1FL' 'LS1HL' 'LM1' 'LM2' 'LCG1' 'LCG2'...
%     'RCG1' 'RCG2' 'RM2' 'RM1' 'RCPU' 'RS1HL' 'RS1FL' 'RS1DZ' 'RS1BF' 'RS2'};
fusroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2'...
    'RCg1' 'RCg2' 'RM2' 'LS1BF' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};
% 
% fmriroi = {'LS2' 'LS1BF' 'LCPU' 'LS1DZ' 'LS1FL' 'LS1HL' 'LM1' 'LM2' 'LCG1' 'LCG2'...
%     'RCG1' 'RCG2' 'RM2' 'RM1' 'RCPU' 'RS1HL' 'RS1FL' 'RS1DZ' 'RS1BF' 'RS2'};
fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
    'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

fmri = load('peak_fMRI.mat'); fmri = fmri.tsnormavepeak;
fus = load('peak_fUS.mat'); fus = fus.tsnormavepeak;
fmrits = load('ts_fMRI.mat');
fusts = load('ts_fUS.mat');

figure(10); clf; figure(11); clf; figure(12); clf; figure(13); clf; figure(14); clf;


C = ['b' 'r' 'g' 'c' 'k' 'y' 'm']; Sem = {'square', 'diamond' ,'o'};
Xtot = []; Ytot = [];
for i = 1:20 % # roi
    
    fusname = fusroi{i};
    fmriidx = find(strcmp(fmriroi,fusname) == 1);
    
    figure(10)
    subplot(5,4,i)
    hold on
    XX = []; YY = [];
    for j = 1:3 % # intensity
        
        X = fmri{j,fmriidx};
        Y = fus{j,i};
        for k = 1:size(X,2)
            scatter(X(k),Y(k),'filled',C(k),Sem{j})
        end
%         scatter(X,Y,50,C(j),'filled')
        XX = [XX X];
        YY = [YY Y];
        
    end
    Xtot = [Xtot mean(XX)]; Ytot = [Ytot mean(YY)]; 
    b = polyfit(XX,YY,1);
    f = polyval(b, XX);
    Bbar = mean(YY);
    SStot = sum((YY - Bbar).^2);
    SSreg = sum((f - Bbar).^2);
    SSres = sum((YY - f).^2);
    R2 = 1 - SSres/SStot;
    
    title([fusname ': R^2 = ' num2str(round(R2*100)/100)]);
    
    figure(11)
    subplot(5,4,i)
    hold on
    stdshade(fusts.tsnorm{1,i},0.25,'b-')
    stdshade(fusts.tsnorm{2,i},0.25,'r-')
    stdshade(fusts.tsnorm{3,i},0.25,'g-')
    title(fusname)
    
    figure(12)
    subplot(5,4,i)
    hold on
    stdshade(fmrits.tsnorm{1,i},0.25,'b-')
    stdshade(fmrits.tsnorm{2,i},0.25,'r-')
    stdshade(fmrits.tsnorm{3,i},0.25,'g-')
    title(fusname)
    
    if strcmp(fusname,'LM1')
        figure(13); subplot(2,2,1); hold on
        stdshade(fusts.tsnorm{1,i},0.25,'b-')
        stdshade(fusts.tsnorm{2,i},0.25,'r-')
        stdshade(fusts.tsnorm{3,i},0.25,'g-')
        set(gca,'ylim',[-25 125],'xlim',[0 240]);
        title(fusname)
        
        figure(14); subplot(2,2,1); hold on
        stdshade(fmrits.tsnorm{1,i},0.25,'b-')
        stdshade(fmrits.tsnorm{2,i},0.25,'r-')
        stdshade(fmrits.tsnorm{3,i},0.25,'g-')
        set(gca,'ylim',[-0.5 1.75],'xlim',[0 240]);
        title(fusname)
        
    elseif strcmp(fusname,'RM1')
        figure(13); subplot(2,2,2); hold on
        stdshade(fusts.tsnorm{1,i},0.25,'b-')
        stdshade(fusts.tsnorm{2,i},0.25,'r-')
        stdshade(fusts.tsnorm{3,i},0.25,'g-')
        set(gca,'ylim',[-25 125],'xlim',[0 240]);
        title(fusname)
        
        figure(14); subplot(2,2,2); hold on
        stdshade(fmrits.tsnorm{1,i},0.25,'b-')
        stdshade(fmrits.tsnorm{2,i},0.25,'r-')
        stdshade(fmrits.tsnorm{3,i},0.25,'g-')
        set(gca,'ylim',[-0.5 1.75],'xlim',[0 240]);
        title(fusname)
        
    elseif strcmp(fusname,'LCPu')
        figure(13); subplot(2,2,3); hold on
        stdshade(fusts.tsnorm{1,i},0.25,'b-')
        stdshade(fusts.tsnorm{2,i},0.25,'r-')
        stdshade(fusts.tsnorm{3,i},0.25,'g-')
        set(gca,'ylim',[-25 125],'xlim',[0 240]);
        title(fusname)
        
        figure(14); subplot(2,2,3); hold on
        stdshade(fmrits.tsnorm{1,i},0.25,'b-')
        stdshade(fmrits.tsnorm{2,i},0.25,'r-')
        stdshade(fmrits.tsnorm{3,i},0.25,'g-')
        set(gca,'ylim',[-0.5 1.75],'xlim',[0 240]);
        title(fusname)
        
    elseif strcmp(fusname,'RCPu')
        figure(13); subplot(2,2,4); hold on
        stdshade(fusts.tsnorm{1,i},0.25,'b-')
        stdshade(fusts.tsnorm{2,i},0.25,'r-')
        stdshade(fusts.tsnorm{3,i},0.25,'g-')
        set(gca,'ylim',[-25 125],'xlim',[0 240]);
        title(fusname)
        
        figure(14); subplot(2,2,4); hold on
        stdshade(fmrits.tsnorm{1,i},0.25,'b-')
        stdshade(fmrits.tsnorm{2,i},0.25,'r-')
        stdshade(fmrits.tsnorm{3,i},0.25,'g-')
        set(gca,'ylim',[-0.5 1.75],'xlim',[0 240]);
        title(fusname)
        
    end
    
    
    
    

end
    
    
%%



