% Main data directory
if isunix
    Dir = '/media/bradley/Seagate Backup Plus Drive/';
    slash = '/';
elseif ispc
    Dir = 'D:\';
    slash = '\';
end
cd(Dir)
% addpath(genpath([Dir 'Preprocessing' slash]));
storage = [Dir 'Data_Processed' slash 'fMRI' slash]; if ~exist(storage,'dir'); mkdir(storage); end
% fmri_fold = [Dir 'fMRI' slash 'ofMRI_pre_CW' slash];
addpath(genpath(storage));

% % Total "pre" ofMRI
% base_fold = {'20191107_111952_BEd_preCW_4362998_R1_D5_1_42';
%     '20191108_083343_BEd_preCW_4346075_N_D5_1_47';
%     '20191108_143147_BEd_preCW_4364124_R1_D1_1_51';
%     '20191108_162123_BEd_preCW_4364124_L1_D1_1_52';
%     '20191111_151047_BEd_preCW_4364122_N_D2_1_53';
%     '20191111_192535_BEd_preCW_4364123_N_D2_1_55';
%     '20191111_220250_BEd_preCW_4364123_R1_D2_2_1_57';
%     '20191116_093138_BEd_preCW_4364121_N_D2_1_58';
%     '20191116_115523_BEd_preCW_4364121_R1_D2_1_59';
%     '20191116_133812_BEd_preCW_4364122_R1_D3_1_60';
%     '20191116_170401_BEd_preCW_4364124_N_D1_1_62';
%     };
% 
% SNR = 1;
% x1 = {'0.1'; [10 11]; [14 15]; [14 15]; [13 14 15]; [16 17]; [14 15]; [9 10] ; [13 14]; [17 18]; [9 10]; [12 13]};
% x2 = {'0.5'; [13 14]; [12 13]; [10 11]; [18]; [11 12]; [16 17]; [5 6]; [17 18]; [10 11]; [14 15]; [16 17]}; 
% x3 = {'1.0'; [17 18]; [10 11]; [12 13]; [11 12]; [9 10]; [12 13]; [7 8]; [11 12]; [15 16]; [12 13]; [18 19]};
% x4 = {'1.5'; [15 16]; [8 9]; [8 9]; [16 17]; [15]; [10 11]; [3 4]; [15 16]; [12 14]; [16 17]; [14 15]};
% stimidx = {x1, x2, x3, x4};

% Total "pre" ofMRI
base_fold = {'20191108_083343_BEd_preCW_4346075_N_D5_1_47';
    '20191108_143147_BEd_preCW_4364124_R1_D1_1_51';
    '20191108_162123_BEd_preCW_4364124_L1_D1_1_52';
    '20191111_151047_BEd_preCW_4364122_N_D2_1_53';
    '20191111_220250_BEd_preCW_4364123_R1_D2_2_1_57';
    '20191116_093138_BEd_preCW_4364121_N_D2_1_58';
    '20191116_133812_BEd_preCW_4364122_R1_D3_1_60';
    };

SNR = 1;
x1 = {'0.1'; [14 15]; [14 15]; [14 15]; [16 17]; [9 10] ; [13 14]; [9 10]};
x2 = {'0.5'; [12 13]; [10 11]; [18]; [11 12]; [5 6]; [17 18]; [14 15]}; 
x3 = {'1.0'; [10 11]; [12 13]; [11 12]; [9 10]; [7 8]; [11 12]; [12 13]};
x4 = {'1.5'; [8 9]; [8 9]; [16 17]; [15]; [3 4]; [15 16]; [16 17]};
stimidx = {x1, x2, x3, x4};

% Total "post" ofMRI
% base_fold = {'20191122_095408_BEd_postCW_4346075_N_1_1';
%     '20191122_123846_BEd_postCW_4364143_N_1_2';
%     '20191122_143826_BEd_postCW_4364124_L1_1_3';
%     '20191124_191348_BEd_postCW_4364122_N_1_4';
%     '20191124_204714_BEd_postCW_4364123_R1_1_5';
%     };
% 
% SNR = 1;
% x1 = {'0.1'; [25 26]; [13 14]; [9 10]; [7 8]; [8 9]};
% x2 = {'0.5'; [27 28 29]; [11 12]; [7 8]; [11 12]; [6 7]}; 
% x3 = {'1.0'; [23 24]; [15 16]; [5 6]; [9 10]; [11 13]};
% stimidx = {x1, x2, x3};

% Total "Control" pre ofMRI
% base_fold = {'20200210_161222_BEd_preCW_4419409_N_1_1_63';
%     '20200211_090028_BEd_preCW_4419410_N_1_1_64';
%     };
% 
% SNR = 1;
% x1 = {'0.1'; [17 20]; [14 15]};
% x2 = {'0.5'; [21 22]; [12 13]};
% x3 = {'1.0'; [24 27]; [10 11]};
% stimidx = {x1, x2, x3};

% spm fmri

idx = 1;
% fMRI_Anat_Template(storage,base_fold,idx)
%% Individual fixed effects
% specify mouse session folder
for i = 1:size(base_fold,1)
    
    % Specify reconstructed functional scans from session
    recon_fold = [storage base_fold{i}(1:8) slash base_fold{i}(17:end) slash];
    func_data = dir(recon_fold);
    func_data(~contains({func_data.name},'f')) = [];
    anat_fold = dir(recon_fold);
    anat_fold(~contains({anat_fold.name},'a')) = [];
    
    if ~isempty(anat_fold); anat_fold = [recon_fold anat_fold(1).name slash]; end
    cd(recon_fold)

    param.preproc = 1;
    param.smooth = 2.5; % um *100
    param.Dummy = 10; %# TRs
    param.order = 4; %GLM gamma
    param.template = 1; % User study brain template

    if isequal(SNR,1) % Average runs for higher SNR

        recon_data = dir(recon_fold);
        recon_data(~contains({recon_data.name},'f')) = [];
        for  j = 1:size(recon_data,1)
            recon_opt(j) = str2double(recon_data(j).name(regexp(recon_data(j).name,'\d')));
        end

        for k = 3%size(stimidx,2) % stim intensities
            
            snr_fold = [recon_fold stimidx{k}{1} slash];
            if ~exist(snr_fold,'dir'); mkdir(snr_fold); end
            func_fold = snr_fold;
            func_file = [stimidx{k}{1} '_EPI.nii'];
            
            if ~exist([snr_fold func_file],'file')
                snridx = stimidx{k}{i+1};
                [aa,bb] = ismember(snridx,recon_opt);
                snr_data = recon_data(bb(aa));
                for l = 1:size(snr_data,1)
                    func_file = [recon_fold snr_data(l).name slash snr_data(l).name(1:end-1) '_EPI.nii'];
                    niitmp = load_nii(func_file);
                    niiimg(:,:,:,:,l) = niitmp.img;
                end
                
                cd(snr_fold)
                snrnii.img = mean(niiimg,5);
                snrnii = make_nii(snrnii.img,[1 1 1]);
                for l = 1:size(snrnii.img,4)
                    tmp(:,:,:,l) = flipud(permute(snrnii.img(:,:,:,l),[2 1 3]));
                end
                snrnii.img = tmp;
                snrnii = make_nii(snrnii.img,[1 1 1]);
                clear tmp
                
                func_file = [stimidx{k}{1} '_EPI.nii'];
                save_nii(snrnii, [func_fold func_file]);
                
            else
                snrnii = load_nii([func_fold func_file]);
                cd(snr_fold)
            end
            
            SNR_func_file = [func_fold 'csnr' func_file];
            if exist(SNR_func_file,'file')
                nii = load_untouch_nii(SNR_func_file);
            else
                fMRI_do_preprocess_snr(func_fold,func_file,snrnii,anat_fold,param);
                nii = load_untouch_nii(SNR_func_file);
            end
            [S{i},T{i},R{i}] = compute_snr(nii.img);
            
%             Func = fMRI_do_preprocess(func_fold,func_file,snrnii,anat_fold,param);
            
        end

    else % Analyze each run individually

        for j = size(func_data,1)
            
            func_fold = [recon_fold func_data(j).name slash];
            func_file = [func_data(j).name(1:end-1) '_EPI.nii'];
            full_file = [func_fold func_file];
            cd(func_fold)
            
            nii = load_nii(full_file);
            for k = 1:size(nii.img,4)
                tmp(:,:,:,k) = flipud(permute(nii.img(:,:,:,k),[2 1 3]));
            end
            nii.img = tmp;
            
            SNR_func_file = [func_fold 'csnr' func_file];
            if exist(SNR_func_file,'file')
                nii = load_untouch_nii(SNR_func_file);
            else
                fMRI_do_preprocess_snr(func_fold,func_file,nii,anat_fold,param);
                nii = load_untouch_nii(SNR_func_file);
            end
            [S{i},T{i},R{i}] = compute_snr(nii.img);
            

%             Func = fMRI_do_preprocess(func_fold,func_file,nii,anat_fold,param);

        end

    end

end
save('D:\Data_Processed\fMRI_SNR.mat','S','T','R')
%% Group level fixed effects

clear all
% Main data directory
if isunix
    Dir = '/media/bradley/Seagate Backup Plus Drive/';
    slash = '/';
elseif ispc
    Dir = 'D:\';
    slash = '\';
end
cd(Dir)
% addpath(genpath([Dir 'Preprocessing' slash]));
addpath(genpath('C:\Users\bedelman\Documents\GitHub\ofUSI\'))
storage = [Dir 'Data_Processed' slash 'fMRI' slash]; if ~exist(storage,'dir'); mkdir(storage); end
% fmri_fold = [Dir 'fMRI' slash 'ofMRI_pre_CW' slash];
addpath(genpath(storage));

base_fold = {'20191108_083343_BEd_preCW_4346075_N_D5_1_47';
    '20191111_151047_BEd_preCW_4364122_N_D2_1_53';
    '20191111_220250_BEd_preCW_4364123_R1_D2_2_1_57';
    '20191108_162123_BEd_preCW_4364124_L1_D1_1_52';
    '20191108_143147_BEd_preCW_4364124_R1_D1_1_51';
    '20191116_093138_BEd_preCW_4364121_N_D2_1_58';
    '20191116_133812_BEd_preCW_4364122_R1_D3_1_60';
    }; param.sliceidx = 10; param.maxt = 40; param.descr = 'pre';
atmp = [storage '20191108' slash 'BEd_preCW_4346075_N_D5_1_47' slash '16a' slash 'ave_anatomy.nii'];
mtmp = [storage '20191108' slash 'BEd_preCW_4346075_N_D5_1_47' slash '16a' slash 'ave_anatomy_mask.nii'];
roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_pre_20191205' slash];
fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
    'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

% base_fold = {'20191122_095408_BEd_postCW_4346075_N_1_1';
%     '20191122_123846_BEd_postCW_4364143_N_1_2';
% %     '20191122_143826_BEd_postCW_4364124_L1_1_3';
%     '20191124_191348_BEd_postCW_4364122_N_1_4';
%     '20191124_204714_BEd_postCW_4364123_R1_1_5';
%     }; param.sliceidx = 10; param.maxt = 40; param.descr = 'single';
% atmp = [storage '20191108' slash 'BEd_preCW_4346075_N_D5_1_47' slash '16a' slash 'ave_anatomy.nii'];
% mtmp = [storage '20191108' slash 'BEd_preCW_4346075_N_D5_1_47' slash '16a' slash 'ave_anatomy_mask.nii'];
% roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_pre_20191205' slash];
% fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
%     'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

% base_fold = {'20191122_095408_BEd_postCW_4346075_N_1_1';
%     '20191122_123846_BEd_postCW_4364143_N_1_2';
%     '20191122_143826_BEd_postCW_4364124_L1_1_3';
%     '20191124_191348_BEd_postCW_4364122_N_1_4';
%     '20191124_204714_BEd_postCW_4364123_R1_1_5';
%     }; param.sliceidx = 10; param.maxt = 20; param.descr = 'post';
% atmp = [storage '20191122' slash 'BEd_postCW_4346075_N_1_1' slash '30a' slash 'prepostave_anatomy.nii'];
% mtmp = [storage '20191122' slash 'BEd_postCW_4346075_N_1_1' slash '30a' slash 'ave_anatomy_mask.nii'];
% roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_post_20200212' slash];
% fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
%     'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};
% 
% base_fold = {'20191108_083343_BEd_preCW_4346075_N_D5_1_47';
%     '20191111_220250_BEd_preCW_4364123_R1_D2_2_1_57'; 
%     '20191111_151047_BEd_preCW_4364122_N_D2_1_53';
%     '20191108_143147_BEd_preCW_4364124_R1_D1_1_51';
%     }; param.sliceidx = 10; param.maxt = 40; param.descr = 'ephys';
% atmp = [storage '20191108' slash 'BEd_preCW_4346075_N_D5_1_47' slash '16a' slash 'ave_anatomy.nii'];
% mtmp = [storage '20191108' slash 'BEd_preCW_4346075_N_D5_1_47' slash '16a' slash 'ave_anatomy_mask.nii'];
% roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_pre_20191205' slash];
% fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
%     'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

% base_fold = {'20200210_161222_BEd_preCW_4419409_N_1_1_63';
%     '20200211_090028_BEd_preCW_4419410_N_1_1_64';
%     }; param.sliceidx = 16; param.maxt = 20; param.descr = 'ctrl';
% atmp = [storage '20200211' slash 'BEd_preCW_4419410_N_1_1_64' slash '16a' slash 'ave_anatomy.nii'];
% mtmp = [storage '20200211' slash 'BEd_preCW_4419410_N_1_1_64' slash '16a' slash 'ave_anatomy_mask.nii'];
% roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_ctrl_20200212' slash];
% fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
%     'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

%%
param.order = 4; %GLM gamma
param.template = 1; % User study brain template; must have template for group level
param.Dummy = 10; %# TRs
param.Pthresh = 0.005;
param.fmriroi = fmriroi;
param.atmp = atmp;
param.mtmp = mtmp;
SNR = 1; % Always stim intensity-based analysis
%%
% Load and compile data from individual animals
fmri_organize_data(storage,base_fold,slash,param)
% Fixed effects analysis for each animal
fmri_individual_fixed_effects(storage,base_fold,slash,param)
% Sig active voxel count (+ flow count)
fmri_active_voxel_count(storage,base_fold,slash,param)
% Fixed effects analysis at the group level
fmri_group_fixed_effects(storage,base_fold,slash,param)
% Plot group activation maps
fmri_plot_group_fixed_effects(storage,base_fold,slash,param)
% Extract time series for roi and vascular components thereof
fmri_timeseries_analysis(storage,base_fold,slash,param)
% Plot time series info
fmri_plot_timeseries_analysis(storage,base_fold,slash,param)
% Plot ROI activation (percent, t-stat)
fmri_plot_ROI_info(storage,base_fold,slash,param)
% Regression Analysis on time series info
fmri_regression_analysis(storage,base_fold,slash,param)

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
ind_cont = ind_cont(:,:,:,:,1:3);
reg_Data = [];
for i = 1:size(base_fold,1)
    for j = 1:size(stim,2)-1
        reg_Data = [reg_Data; reshape(ind_cont(:,:,sliceidx,i,j),1,size(ind_cont,1)*size(ind_cont,2))];
    end
end

reg_c = [1;zeros(size(base_fold,1),1)]; % Regression contrast
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

atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy.nii';
anii = load_nii(atmp); anii.img = anii.img(:,:,sliceidx);
mtmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy_mask.nii';
mnii = load_nii(mtmp);
mask = mnii.img(:,:,sliceidx);
% Resize time-series param to anat sampling
reg_p2tail_tmp = imresize(reg_c_p2tail,[size(anii.img,1) size(anii.img,2)]); reg_p2tail_tmp(mask == 0) = 1.1;
reg_t_stat_tmp = imresize(reg_t_stat,[size(anii.img,1) size(anii.img,2)]); reg_t_stat_tmp(mask == 0) = 0;

% rescale dynamic range of anat for vis
Pthresh = 0.05;
anat = anii.img; anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
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
    for j = 1:size(stim,2)-1
        anov_DATA = [anov_DATA; reshape(ind_cont(:,:,sliceidx,i,j),1,size(ind_cont,1)*size(ind_cont,2))];
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

atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy.nii';
anii = load_nii(atmp); anii.img = anii.img(:,:,sliceidx);
mtmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy_mask.nii';
mnii = load_nii(mtmp);
mask = mnii.img(:,:,sliceidx);
% Resize time-series param to anat sampling
anov_F_tmp = imresize(anov_F,[size(anii.img,1) size(anii.img,2)]); anov_F_tmp(mask == 0) = plotmin - maxf;
anov_p_tmp = imresize(anov_p,[size(anii.img,1) size(anii.img,2)]); anov_p_tmp(mask == 0) = 1.1;
ttest1_p_tmp = imresize(ttest1_p,[size(anii.img,1) size(anii.img,2)]); ttest1_p_tmp(mask == 0) = 1.1;
ttest1_t_stat_tmp = imresize(ttest1_t_stat,[size(anii.img,1) size(anii.img,2)]); ttest1_t_stat_tmp(mask == 0) = plotmin - maxt;
ttest2_p_tmp = imresize(ttest2_p,[size(anii.img,1) size(anii.img,2)]); ttest2_p_tmp(mask == 0) = 1.1;
ttest2_t_stat_tmp = imresize(ttest2_t_stat,[size(anii.img,1) size(anii.img,2)]); ttest2_t_stat_tmp(mask == 0) = plotmin - maxt;
ttest3_p_tmp = imresize(ttest3_p,[size(anii.img,1) size(anii.img,2)]); ttest3_p_tmp(mask == 0) = 1.1;
ttest3_t_stat_tmp = imresize(ttest3_t_stat,[size(anii.img,1) size(anii.img,2)]); ttest3_t_stat_tmp(mask == 0) = plotmin - maxt;

% rescale dynamic range of anat for vis
Pthresh = 0.05;
anat = anii.img; anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
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









