% Main data directory
if isunix
    Dir = '/media/bradley/Seagate Backup Plus Drive/';
    slash = '/';
elseif ispc
    Dir = 'D:\';
    slash = '\';
end
cd(Dir)
addpath(genpath([Dir 'Preprocessing' slash]));
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

spm fmri

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

        for k = 1:3%size(stimidx,2) % stim intensities
            
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
            
            Func = fMRI_do_preprocess(func_fold,func_file,snrnii,anat_fold,param);
            
        end

    else % Analyze each run individually

        for j = 1:size(func_data,1)

            func_fold = [recon_fold func_data(j).name slash];
            func_file = [func_data(j).name(1:end-1) '_EPI.nii'];
            full_file = [func_fold func_file];
            cd(func_fold)

            nii = load_nii(full_file);
            for k = 1:size(nii.img,4)
                tmp(:,:,:,k) = flipud(permute(nii.img(:,:,:,k),[2 1 3]));
            end
            nii.img = tmp;

            Func = fMRI_do_preprocess(func_fold,func_file,nii,anat_fold,param);

        end

    end

end
% cd(fmri_fold)
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
addpath(genpath([Dir 'Preprocessing' slash]));
storage = [Dir 'Data_Processed' slash 'fMRI' slash]; if ~exist(storage,'dir'); mkdir(storage); end
% fmri_fold = [Dir 'fMRI' slash 'ofMRI_pre_CW' slash];
addpath(genpath(storage));

param.order = 4; %GLM gamma
param.template = 1; % User study brain template; must have template for group level
SNR = 1; % Always stim intensity-based analysis

base_fold = {'20191108_083343_BEd_preCW_4346075_N_D5_1_47';
    '20191111_151047_BEd_preCW_4364122_N_D2_1_53';
    '20191111_220250_BEd_preCW_4364123_R1_D2_2_1_57';
    '20191108_162123_BEd_preCW_4364124_L1_D1_1_52';
    '20191108_143147_BEd_preCW_4364124_R1_D1_1_51';
    '20191116_093138_BEd_preCW_4364121_N_D2_1_58';
    '20191116_133812_BEd_preCW_4364122_R1_D3_1_60';
    }; sliceidx = 10; maxt = 40; descr = 'pre';
atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy.nii';
mtmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy_mask.nii';
roifolder = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/segment_ofUS/fmri_pre_20191205/';
fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
    'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

% base_fold = {'20191122_095408_BEd_postCW_4346075_N_1_1';
%     '20191122_123846_BEd_postCW_4364143_N_1_2';
%     '20191122_143826_BEd_postCW_4364124_L1_1_3';
%     '20191124_191348_BEd_postCW_4364122_N_1_4';
%     '20191124_204714_BEd_postCW_4364123_R1_1_5';
%     }; sliceidx = 10; maxt = 20; descr = 'post';
% atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191122/BEd_postCW_4346075_N_1_1/30a/ave_anatomy.nii';
% mtmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191122/BEd_postCW_4346075_N_1_1/30a/ave_anatomy_mask.nii';
% roifolder = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/segment_ofUS/fmri_post_20200212/';
% fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
%     'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

% base_fold = {'20191108_083343_BEd_preCW_4346075_N_D5_1_47';
%     '20191111_220250_BEd_preCW_4364123_R1_D2_2_1_57'; 
%     '20191111_151047_BEd_preCW_4364122_N_D2_1_53';
%     '20191108_143147_BEd_preCW_4364124_R1_D1_1_51';
%     }; sliceidx = 10; maxt = 40; descr = 'ephys';
% atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy.nii';
% mtmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy_mask.nii';
% roifolder = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/segment_ofUS/fmri_pre_20191205/';
% fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
%     'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

% base_fold = {'20200210_161222_BEd_preCW_4419409_N_1_1_63';
%     '20200211_090028_BEd_preCW_4419410_N_1_1_64';
%     }; sliceidx = 16; maxt = 20; descr = 'ctrl';
% atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20200211/BEd_preCW_4419410_N_1_1_64/16a/ave_anatomy.nii';
% mtmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20200211/BEd_preCW_4419410_N_1_1_64/16a/ave_anatomy_mask.nii';
% roifolder = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/segment_ofUS/fmri_ctrl_20200212/';    
% fmriroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2' 'LS1BF'...
%     'RCg1' 'RCg2' 'RM2' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};


stim = {'0.1', '0.5', '1.0'};

ts = cell(3,1); tsnorm = cell(3,1); pts = cell(3,1);
tsnormave = cell(3,1); tsnormavepeak = cell(3,1); tsnormaveAUC = cell(3,1);
ind_cont = [];
ind_pts = cell(1,size(base_fold,1),3);
ind_p2_roi = cell(1,size(base_fold,1),3);
ind_p2_roi_prop = cell(1,size(base_fold,1),3);
ind_t_stat = zeros(35,80,size(base_fold,1),3);
ind_p2tail = zeros(35,80,size(base_fold,1),3);

Pthresh = 0.005;

for  ii = 1:size(stim,2)
    ii
    ind_data = cell(0); x = cell(0); filt = cell(0); D1 = []; D2 = []; D3 = []; grp_data = [];
    for j = 1:size(base_fold,1)
    
        % Specify reconstructed functional scans from session
        recon_fold = [storage base_fold{j}(1:8) slash base_fold{j}(17:end) slash];
        func_data = dir(recon_fold);
        func_data(~contains({func_data.name},stim{ii})) = [];
        func_file = [func_data.folder slash func_data.name slash 'template' slash 'Func.mat'];

        if exist(func_file,'file')

            load(func_file)

            ind_data{j} = FUNC.Data.smooth.img;
            x{j} = FUNC.GLM.X(:,1:5);
            filt{j} = FUNC.GLM.X(:,6:12);
            D1 = blkdiag(D1,x{j}(:,1:4));
            D2 = blkdiag(D2,x{j}(:,5));
            D3 = blkdiag(D3,filt{j});
            grp_data = cat(4,grp_data,FUNC.Data.smooth.img);
    
        end
        
    end
    
    % Define t-ditributions
    tdist2T = @(t,DF) (1-betainc(DF/(DF+t^2),DF/2,0.5)); % 2-tailed

    % Individual animal fixed effect analysis
    ind_c = [1 0 0 0 0 zeros(1,7)]';
    ind_beta = zeros(size(ind_data{1},1),size(ind_data{1},2),size(ind_data{1},3),size(ind_c,1),size(ind_data,2));
%     ind_t_stat = zeros(size(ind_data{1},1),size(ind_data{1},2),size(ind_data{1},3),size(ind_data,2));
%     ind_p2tail = zeros(size(ind_data{1},1),size(ind_data{1},2),size(ind_data{1},3),size(ind_data,2));
    for j = 1:size(base_fold,1)
        
        DF = size(ind_data{j},3) - 2;
        ind_D = [x{j} filt{j}];
        
        for k = 1:size(ind_data{j},1)
            for l = 1:size(ind_data{j},2)
                for m = 1:size(ind_data{j},3)
                    Y = squeeze(ind_data{j}(k,l,m,:))';

                    ind_beta_hat = inv(ind_D'*ind_D)*ind_D'*Y';
                    ind_beta(k,l,m,:,j) = reshape(ind_beta_hat,1,1,size(ind_c,1));
                    ind_Var_e = (Y'-ind_D*ind_beta_hat)'*(Y'-ind_D*ind_beta_hat)/DF;

                    %Hypothesis testing; Compute the t statistic
                    ind_t_stat(k,l,m,j,ii) = ind_c'*ind_beta_hat/sqrt(ind_Var_e*ind_c'*inv(ind_D'*ind_D)*ind_c);

                    if isnan(ind_t_stat(k,l,m,j,ii))
                        ind_p2tail(k,l,m,j,ii) = 1;
                    else
                        ind_p2tail(k,l,m,j,ii) = 1 - tdist2T(ind_t_stat(k,l,m,j,ii),DF);
                    end
                end
            end
        end
        
        for k = 1:size(ind_beta,3)
            
            ind_cont(:,:,k,j,ii) = sum(repmat(reshape(ind_c,1,[],size(ind_c,1)),...
                size(ind_data{j},1),size(ind_data{j},2)).*squeeze(ind_beta(:,:,k,:,j)),3);
        
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    for i = 1:size(base_fold,1)
        ttt = ind_p2tail(:,:,sliceidx,i,ii);
        subplot(2,size(base_fold,1),i); imagesc(ttt); caxis([0 1])
        ttt(ttt > 0.05) = 1;
        subplot(2,size(base_fold,1),size(base_fold,1)+i); imagesc(ttt);
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract number of significantly active voxels per ROI (only for slice of interest)
    roi = dir(roifolder);
    roi(contains({roi.name},'.fig')) = [];
    roi(~contains({roi.name},'R')) = [];
    
    for i = 1:size(roi,1)
        roifile = [roifolder roi(i).name];
        load(roifile)
        cdata = imresize(ROI.cdata,.2);
        for j = 1:size(base_fold,1)
            
            ind_pts{i,j,ii} = find(cdata); % Total points in roi
            
            tmp = ind_t_stat(:,:,sliceidx,j,ii);
            ind_t_stat_ave{i,j,ii} = mean(tmp(ind_pts{i,j,ii}));
            
            tmp = ind_p2tail(:,:,sliceidx,j,ii);
            tmp = tmp(ind_pts{i,j,ii});
            ind_p2_roi{i,j,ii} = sum(tmp < 0.05); % # of sig points
            ind_p2_roi_prop{i,j,ii} = ind_p2_roi{i,j,ii}/size(ind_pts{i,j,ii},1); % proportion of ROI sig activated
            
        end
    end
    
    % Group fixed effect analysis
    grp_D = [D1 D2 D3];
    nTRs = size(grp_D,1);
    grp_c = [repmat([1 0 0 0],[1 size(base_fold,1)]), zeros(1,size(base_fold,1)), zeros(1,7*size(base_fold,1))]';
    grp_beta = zeros(size(grp_data,1),size(grp_data,2),size(grp_data,3),size(grp_c,1));
    grp_t_stat = zeros(size(grp_data,1),size(grp_data,2),size(grp_data,3));
    grp_p2tail = zeros(size(grp_data,1),size(grp_data,2),size(grp_data,3));
    for i = 1:size(grp_data,1)
        for j = 1:size(grp_data,2)
            for k = 1:size(grp_data,3)

                Y = squeeze(grp_data(i,j,k,:))';

                DF = nTRs - 2;

                beta_hat=inv(grp_D'*grp_D)*grp_D'*Y';
                grp_beta(i,j,k,:) = reshape(beta_hat,1,1,size(grp_c,1));
                Var_e=(Y'-grp_D*beta_hat)'*(Y'-grp_D*beta_hat)/DF;

                %Hypothesis testing; Compute the t statistic
                grp_t_stat(i,j,k)=grp_c'*beta_hat/sqrt(Var_e*grp_c'*inv(grp_D'*grp_D)*grp_c);

                if isnan(grp_t_stat(i,j,k))
                    grp_p2tail(i,j,k) = 1;
                else
                    grp_p2tail(i,j,k) = 1 - tdist2T(grp_t_stat(i,j,k),DF);
                end

            end
        end
    end
    
    % Apply correction
    [grp_c_p2tail, c_alpha, h] = fwer_sidak(reshape(grp_p2tail,size(grp_p2tail,1)*size(grp_p2tail,2)*size(grp_p2tail,3), []), 0.05);
    grp_c_p2tail = reshape(grp_c_p2tail,size(grp_p2tail,1),size(grp_p2tail,2),size(grp_p2tail,3));
    [h, crit_p, adj_ci_cvrg, grp_c_p2tail]=fdr_bh(grp_p2tail,Pthresh,'pdep','no');
%     c_p2tail = p2tail;
    
    h1 = figure(10); subplot(2,1,1); 
    b=bar(.5:1:size(grp_c,1),grp_c,'k','barwidth',.5); grid minor; title('Design Matrix')
    set(gca,'position',[.13 .80 .775 .15],'ylim',[-2 2],'xlim',[0 size(grp_c,1)],...
        'xtick',0:1:12,'xticklabel',[],'yticklabel',[])
    subplot(2,1,2); image(grp_D*64);
    set(gca,'position',[.13 .1 .775 .65],'xticklabel',[],'yticklabel',[1 size(grp_D,1)],'ytick',[1 size(grp_D,2)])
    ylabel('Acq Number');
    colormap gray

    plotmin = -100; plotmax = 0;
    Cbar = [gray(abs(plotmax - plotmin)*100); fireice(2*maxt*100)]; CAX = [plotmin - maxt maxt];
    Cbarp = [jet(round(1.1*100));[0 0 0]];
    h2 = figure(12);
    
    anii = load_nii(atmp);
    mnii = load_nii(mtmp);
    mask = mnii.img;
    % mask = ones(size(anii.img));
    % Resize time-series param to anat sampling
    grp_c_p2tail_tmp = imresize3(grp_c_p2tail,[size(anii.img,1) size(anii.img,2) size(anii.img,3)]);
    grp_c_p2tail_tmp(mask == 0) = 1.1;
    grp_t_stat_tmp = imresize3(grp_t_stat,[size(anii.img,1) size(anii.img,2) size(anii.img,3)]);
    grp_t_stat_tmp(mask == 0) = 0;

    % rescale dynamic range of anat for vis
    anat = anii.img;
    anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
    anat = anat - maxt;
    grp_anat = anat; grp_anat(grp_c_p2tail_tmp < Pthresh) = grp_t_stat_tmp(grp_c_p2tail_tmp < Pthresh); grp_anat(mask == 0) = 0;
    for i = 1:size(grp_anat,3)
        figure(100); subplot(6,4,i);
        imagesc(grp_c_p2tail_tmp(:,:,i)); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbarp); caxis([0 1.1])

        figure(101); subplot(6,4,i);
        imagesc(grp_anat(:,:,i)); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbar); caxis(CAX); 
    end
    
    figure(h2);
    subplot(4,2,1+2*(ii-1)); imagesc(grp_c_p2tail_tmp(:,:,sliceidx));
    set(gca,'xticklabel',[],'yticklabel',[]); colormap(gca,Cbarp); caxis([0 1.1]); title(['p-val 2tail - ' stim{ii} 'mW'])
    subplot(4,2,2+2*(ii-1)); imagesc(grp_anat(:,:,sliceidx));
    set(gca,'xticklabel',[],'yticklabel',[]); colormap(gca,Cbar); caxis(CAX); title(['t-val 2tail - ' stim{ii} 'mW'])
    
    roifolder = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/tmp_segment_ofUS/fmri_20191205/';
    roi = dir(roifolder);
    roi(contains({roi.name},'.fig')) = [];
    roi(~contains({roi.name},'R')) = [];
    
    % Time series analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROI
    start = [30 70 110 150 190];
    map = zeros(size(imresize(anat(:,:,sliceidx),[35,80])));
    for i = 1:size(roi,1)
        roifile = [roifolder roi(i).name];
        load(roifile)
        cdata = imresize(ROI.cdata,.2);
        pts {ii,i}= find(cdata);
        map(pts{ii,i}) = i;
        
        for j = 1:size(base_fold,1)
            datatmp = reshape(squeeze(ind_data{j}(:,:,sliceidx,:)),[size(ind_data{j},1)*size(ind_data{j},2),size(ind_data{j},4)]);
            [b,a] = butter(4,[0.0005 0.1]/((1/1.5)/2));
            
            ts{ii,i}(j,:) = sum(datatmp(pts{ii,i},:));
            meanbase = mean(ts{ii,i}(j,1:20));
            tsnorm{ii,i}(j,:) = (ts{ii,i}(j,:) - meanbase)/meanbase*100;
            tsnormf{ii,i}(j,:) = filtfilt(b,a,tsnorm{ii,i}(j,:));

            ts_normpeaks = zeros(size(start,2),size(-10:28,2));
            for k = 1:size(start,2)
                ts_normpeaks(k,:) = tsnormf{ii,i}(j,start(k)-10:start(k)+28);
            end
            tmp = ts_normpeaks'-repmat(mean(ts_normpeaks',1),[39 1]);
            tmp = tmp - repmat(mean(tmp(1:10,:),1),[39,1]);
            figure(1000); subplot(2,7,j); plot(tmp)
            subplot(2,7,j+7); plot(mean(tmp,2))
            tsnormave{ii,i}(j,:) = mean(tmp,2)';
            tsnormavepeak{ii,i}(j) = max(mean(tmp'));
            tsnormaveAUC{ii,i}(j) = sum(mean(tmp'));
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIXEL
    start = [30 70 110 150 190];
    [b,a] = butter(4,[0.0005 0.1]/((1/1.5)/2));
    for i = 1:size(base_fold,1)
        datatmp = squeeze(ind_data{i}(:,:,sliceidx,:));
        
        for j = 1:size(ind_data{1},1)
            for k = 1:size(ind_data{1},2)
                meanbase_pix = mean(datatmp(j,k,:));
                tsnorm_pix(j,k,:,i,ii) = (datatmp(j,k,:) - meanbase_pix)/meanbase_pix*100;
                tsnorm_pix(j,k,:,i,ii) = filtfilt(b,a,tsnorm_pix(j,k,:,i,ii));
                
                ts_normpeaks = zeros(size(start,2),size(-10:28,2));
                for l = 1:size(start,2)
                    ts_normpeaks(l,:) = squeeze(tsnorm_pix(j,k,start(l)-10:start(l)+28,i,ii));
                end
                tmp = ts_normpeaks'-repmat(mean(ts_normpeaks',1),[39 1]);
                tmp = tmp - repmat(mean(tmp(1:10,:),1),[39,1]);
                tsnormave_pix(j,k,:,i,ii) = mean(tmp,2)';
                tsnormavepeak_pix(j,k,i,ii) = max(mean(tmp'));
                tsnormaveAUC_pix(j,k,i,ii) = sum(mean(tmp'));
            end
        end
    end
    
    maxp = 1.75;
    anattmp = -2*maxp*ones(size(imresize(imresize(anat(:,:,sliceidx),[35 80]),5)));
    for i = 1:size(roi,1)
        roifile = [roifolder roi(i).name];
        load(roifile)
        pts2{i}= find(ROI.cdata);
        anattmp(pts2{i}) = mean(tsnormavepeak{ii,i},2); mean(tsnormavepeak{ii,i},2)
    end
%     maxp = max(max(anattmp)); maxp = 1.5;
    anatroi = imresize(imresize(anat(:,:,sliceidx) + maxt,[35 80]),5);
    anatroi = anatroi - maxp;
    anatroi(anattmp ~= -2*maxp) = anattmp(anattmp ~= -2*maxp);
    
    Cbar = [gray(abs(plotmax - plotmin)*100); jet(2*maxp*100)]; CAX = [plotmin - maxp maxp];
    figure; subplot(1,2,1); imagesc(anatroi); caxis(CAX);
    subplot(1,2,2); imagesc(anattmp); colormap(Cbar); caxis(CAX);
    
end

% Save a few variables: time series peak, AUC
save([storage 'fMRI_time_series_peak_' descr '.mat'],'tsnormavepeak','fmriroi');
save([storage 'fMRI_time_series_AUC_' descr '.mat'],'tsnormaveAUC','fmriroi');

figure(500); clf
for j = 1:20
    subplot(5,4,j); hold on
    
    tmp1 = tsnormave{1,j};
    stdshade(tmp1,0.25,[227 140 187]/255)
    tmp2 = tsnormave{2,j};
    stdshade(tmp2,0.25,[176 41 140]/255)
    tmp3 = tsnormave{3,j};
    stdshade(tmp3,0.25,[86 30 88]/255)

    set(gca,'ylim',[-.5 2])
    y1=get(gca,'ylim');
    plot([10 10],y1,'k'); plot([18 18],y1,'k')
    title(fmriroi{j})
end

% Proportion ROI activated
Prop = cell2mat(ind_p2_roi_prop(:,:,1:3));
PropM = squeeze(mean(Prop,2));
PropSem = squeeze(std(Prop,0,2)/sqrt(size(Prop,2)));

figure;
group = cell(size(base_fold,1),1); group(:) = {'1'};
for i = 1:size(roi,1)
    for j = 1:size(base_fold,1)
        Prop2(j,:,i) = squeeze(Prop(i,j,:))';
    end
    
    subplot(5,4,i); hold on
    parallelcoords(Prop2(:,:,i),'Group',group,'Labels',stim(1:3),'color','k');
    errorbar(1:3,mean(Prop2(:,:,i),1),std(Prop2(:,:,i),1)/sqrt(size(Prop2(:,:,i),1)),...
        '-o','markersize',3,'linewidth',1.5,'color','r')
    title(fmriroi{i})

    legend off; set(gca,'ylabel',[],'ylim',[0 1.25],'xlim',[0.75 3.25])
    
    [H,P12] = ttest(Prop2(:,1,i),Prop2(:,2,i)); text(1.35,1,num2str(P12,'%1.2f'))
    [H,P23] = ttest(Prop2(:,2,i),Prop2(:,3,i)); text(2.35,1,num2str(P23,'%1.2f'))
    [H,P13] = ttest(Prop2(:,1,i),Prop2(:,3,i)); text(1.85,1.25,num2str(P13,'%1.2f'))
end
suptitle('Proportion ROI Activated')

ctrs = 1:20;
data = PropM;
figure
hBar{1} = bar(ctrs, data);
hold on
for k1 = 1:size(PropM,2)
    ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
    ydt(k1,:) = hBar{1}(k1).YData;
end
errorbar(ctr, ydt, zeros(3,20), PropSem', zeros(3,20), zeros(3,20), '.k',...
    'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
hBar{2} = bar(ctrs, data); delete(hBar{1})
set(gca,'ylim',[0 1],'xtick',[])
hold off
suptitle('Proportion ROI Activated')

% Ave t-stat
Tval = cell2mat(ind_t_stat_ave);
TvalM = squeeze(nanmean(Tval,2)); TvalM = TvalM(:,1:3);
TvalSem = squeeze(nanstd(Tval,0,2)/sqrt(size(Tval,2))); TvalSem = TvalSem(:,1:3);
ctrs = 1:20;
data = TvalM;
figure; hBar{1} = bar(ctrs, data); hold on
clear ctr ydt
for k1 = 1:size(TvalM,2)
    ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
    ydt(k1,:) = hBar{1}(k1).YData;
end
errorbar(ctr, ydt, zeros(3,20), TvalSem', zeros(3,20), zeros(3,20),'.k',...
    'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
hBar{2} = bar(ctrs, data); delete(hBar{1})
set(gca,'ylim',[-.5 3],'xtick',1:20,'xticklabel',fmriroi)
hold off; suptitle('Average T-value')

clear StatROIt
StatROIt(1,1:3) = {'subjid','intensity','pairInt'};
id = repmat(1:size(base_fold,1),3,1); StatROIt(2:2+3*size(base_fold,1)-1,1) = num2cell(id(:));
StatROIt(2:2+3*size(base_fold,1)-1,2) = num2cell(repmat([1;2;3],size(base_fold,1),1));
StatROIt(2:2+3*size(base_fold,1)-1,3) = num2cell(repmat([1;2;3],size(base_fold,1),1));
for  i = 1:20
        tmp1 = cell2mat(squeeze(ind_t_stat_ave(i,:,:)))';
        StatROIt(1,3+i) = fmriroi(i);
        StatROIt(2:2+3*size(base_fold,1)-1,3+i) = num2cell(tmp1(:));
end
        
% Write stat file for R - ROI average T values
StatROIt(:,[9,23]) = [];
fid = fopen(['/home/bradley/Dropbox/ROI_Tval_fMRI_' descr '.txt'],'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatROIt{1,:});
for K = 2:size(StatROIt,1)
    fprintf(fid, '%.0f %.2f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatROIt{K,:});
end
fclose(fid)




%%%%%%%% REGRESSION
BB = zeros(size(ind_t_stat,1),size(ind_t_stat,2));
RR = zeros(size(ind_t_stat,1),size(ind_t_stat,2));
PP = zeros(size(ind_t_stat,1),size(ind_t_stat,2));
for i = 1:size(ind_t_stat,1)
    for j = 1:size(ind_t_stat,2)
        
        tmp1 = squeeze(ind_t_stat(i,j,sliceidx,:,1));
        tmp2 = squeeze(ind_t_stat(i,j,sliceidx,:,2));
        tmp3 = squeeze(ind_t_stat(i,j,sliceidx,:,3));
        
        tmp4 = 0.1*ones(size(base_fold,1),1);
        tmp5 = 0.5*ones(size(base_fold,1),1);
        tmp6 = 1.0*ones(size(base_fold,1),1);

        T1 = [tmp1;tmp2;tmp3];
        T2 = [tmp4;tmp5;tmp6];
        T2(isnan(T1)) = [];
        T1(isnan(T1)) = [];
        
        [B,BINT,R,RINT,STATS] = regress(T1,[ones(size(T2,1),1),T2]);
        
        BB(i,j) = B(2);
        RR(i,j) = STATS(1);
        PP(i,j) = STATS(3);
        
        if i == 14 && j == 44
            
            figure; scatter(T2,T1); hold on
            set(gca,'xlim',[-.5 1.5],'ylim',[-2 5]);
            plot(reshape(T2,[7,3])',reshape(T1,[7,3])')
            scatter([0.1 0.5 1.0], mean(reshape(T1,[7,3]),1),50,'k','filled')
            plot([0.1 0.5 1.0],mean(reshape(T1,[7,3]),1),'k','linewidth',2)
            
        end
        
    end
end


% rescale dynamic range of anat for vis
maxR = max(RR(:)); maxB = max(abs(BB(:))); maxB = 15; plotmin = -100; plotmax = 0;
CbarR = [gray(abs(plotmax - plotmin)*100); fireice(2*maxR*100)]; CAXR = [plotmin - maxR maxR];
CbarB = [gray(abs(plotmax - plotmin)*100); fireice(2*maxB*100)]; CAXB = [plotmin - maxB maxB];

BB2 = imresize(BB,[256,256]); BB2(mask(:,:,sliceidx) == 0) = CAXB(1);
RR2 = imresize(RR,[256,256]); RR2(mask(:,:,sliceidx) == 0) = CAXB(1);
PP2 = imresize(PP,[256,256]);
anat = anii.img(:,:,sliceidx);
anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin; 
anatRR = anat - maxR; anatBB = anat - maxB;
anatBB(PP2 < .05) = BB2(PP2 < .05); anatBB(mask(:,:,sliceidx) == 0) = CAXB(1);
anatRR(PP2 < .05) = RR2(PP2 < .05); anatRR(mask(:,:,sliceidx) == 0) = CAXR(1);

figure;
subplot(2,2,1); imagesc(BB2); colormap(gca,CbarB); caxis(CAXB); title('fMRI - Regression Beta');
subplot(2,2,2); imagesc(RR2); colormap(gca,CbarR); caxis(CAXR); title('fMRI - Regression Corr'); 
BB2(PP2>0.05) = 0; RR2(PP2>0.05) = 0;
subplot(2,2,3); imagesc(anatBB); colormap(gca,CbarB); caxis(CAXB); title('fMRI - Regression Beta Thresh');
subplot(2,2,4); imagesc(anatRR); colormap(gca,CbarR); caxis(CAXR); title('fMRI - Regression Corr Thresh');



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









