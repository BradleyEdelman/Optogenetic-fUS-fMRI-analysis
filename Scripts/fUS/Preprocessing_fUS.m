%% Bulk preprocessing script

% Main data directory
if isunix
    Dir = '/media/bradley/Seagate Backup Plus Drive/';
    slash = '/';
elseif ispc
    Dir = 'D:\';
    slash = '\';
end
addpath(genpath([Dir 'Preprocessing' slash]));
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
        
        for j = 1:size(func_data,1)
            
            func_fold = [recon_fold func_data(j).name(end-16:end-9) slash];
            if ~exist(func_fold,'dir'); mkdir(func_fold); end
            cd(func_fold)
            
            anat_fold = recon_fold;
            
            FUNC = fUS_do_preprocess(func_fold,func_data(j),anat_fold,param);
                    
        end
    end
%     pause
end


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
addpath(genpath([Dir 'Preprocessing' slash]));
storage = [Dir 'Data_Processed' slash 'fUS' slash ]; if ~exist(storage,'dir'); mkdir(storage); end
addpath(genpath(storage));

param.order = 4; %GLM gamma
param.template = 1; % User study brain template; must have template for group level
param.Dummy = 10; %# TRs
SNR = 1; % Always stim intensity-based analysis

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

stim = {'0.1', '0.5', '1.0'};

ts = cell(3,1); tsnorm = cell(3,1); pts = cell(3,1);
tsnormave = cell(3,1); tsnormavepeak = cell(3,1);
ind_cont = [];
ind_pts = cell(1,size(base_fold,1),3);
ind_p2_roi = cell(1,size(base_fold,1),3);
ind_p2_roi_prop = cell(1,size(base_fold,1),3);
grp_ave_t_val = cell(1,3);
ind_t_stat = zeros(110,80,size(base_fold,1),3);
ind_p2tail = zeros(110,80,size(base_fold,1),3);

Pthresh = 0.005; maxt = 40;

fusroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2'...
    'RCg1' 'RCg2' 'RM2' 'LS1BF' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};

for  ii = 1:size(stim,2)
    ii
    
    % Load individual animal data and construct group-level design matrix
    ind_data = cell(0); x = cell(0); filt = cell(0); D1 = []; D2 = []; D3 = []; grp_data = [];
    for j = 1:size(base_fold,1)
        
        % Specify reconstructed functional scans from session
        recon_fold = [storage base_fold{j}(1:8) slash base_fold{j}(10:end) slash];
        snr_fold = [recon_fold stim{ii} slash];
        func_file = [snr_fold 'template' slash stim{ii} '_Func.mat'];
        
        if exist(func_file,'file')

            load(func_file)
            
            ind_data{j} = imresize3(FUNC.Data.pdi(:,:,param.Dummy+1:end), [110 80 230]);
            x{j} = FUNC.GLM.X(:,1:5);
            filt{j} = FUNC.GLM.X(:,6:12);
            D1 = blkdiag(D1,x{j}(:,1:4));
            D2 = blkdiag(D2,x{j}(:,5));
            D3 = blkdiag(D3,filt{j});
            grp_data = cat(4,grp_data,imresize3(FUNC.Data.pdi(:,:,param.Dummy+1:end), [110 80 230]));
            
        end
        
    end
    
    % Define t-ditributions
    tdist2T = @(t,DF) (1-betainc(DF/(DF+t^2),DF/2,0.5)); % 2-tailed
    tdist1T = @(t,DF) 1-(1-tdist2T(t,DF))/2; % 1-tailed
    
    % Individual animal fixed effect analysis
    ind_c = [1 0 0 0 0 zeros(1,7)]';
    ind_beta = zeros(size(ind_data{1},1),size(ind_data{1},2),size(ind_c,1),size(ind_data,3));
%     ind_t_stat = zeros(size(ind_data{1},1),size(ind_data{1},2),size(ind_data,3));
%     ind_p2tail = zeros(size(ind_data{1},1),size(ind_data{1},2),size(ind_data,3));
    for j = 1:size(base_fold,1)
        DF = size(ind_data{j},3) - 2;
        ind_D = [x{j} filt{j}];
        for k = 1:size(ind_data{j},1)
            for l = 1:size(ind_data{j},2)
                
                Y = squeeze(ind_data{j}(k,l,:))';
                
                ind_beta_hat = inv(ind_D'*ind_D)*ind_D'*Y';
                ind_beta(k,l,:,j) = reshape(ind_beta_hat,1,1,size(ind_c,1));
                ind_Var_e = (Y'-ind_D*ind_beta_hat)'*(Y'-ind_D*ind_beta_hat)/DF;
                
                %Hypothesis testing; Compute the t statistic
                ind_t_stat(k,l,j,ii)=ind_c'*ind_beta_hat/sqrt(ind_Var_e*ind_c'*inv(ind_D'*ind_D)*ind_c);
                
                if isnan(ind_t_stat(k,l,j,ii))
                    ind_p2tail(k,l,j,ii) = 1;
                else
                    ind_p2tail(k,l,j,ii) = 1 - tdist2T(ind_t_stat(k,l,j,ii),DF);
                end
    
            end
        end
        
        ind_cont(:,:,j,ii) = sum((repmat(reshape(ind_c,1,[],size(ind_c,1)),...
            size(ind_data{j},1),size(ind_data{j},2))).*ind_beta(:,:,:,j),3);
        
    end
    
    % Extract number of significantly active voxels per ROI
    roifolder = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/tmp_segment_ofUS/fus_20191205/';
    roi = dir(roifolder);
    roi(contains({roi.name},'.fig')) = [];
    roi(~contains({roi.name},'R')) = [];
    for i = 1:size(roi,1)
        roifile = [roifolder roi(i).name];
        load(roifile)
        cdata = imresize(ROI.cdata,.5);
        for j = 1:size(base_fold,1)
            ind_pts{i,j,ii} = find(cdata); % Total points in roi
            
            tmp = ind_t_stat(:,:,j,ii);
            ind_t_stat_ave{i,j,ii} = mean(tmp(ind_pts{i,j,ii}));
            
            tmp = ind_p2tail(:,:,j,ii);
            tmp = tmp(ind_pts{i,j,ii});
            ind_p2_roi{i,j,ii} = sum(tmp < 0.05); % # of sig points
            ind_p2_roi_prop{i,j,ii} = ind_p2_roi{i,j,ii}/size(ind_pts{i,j,ii},1); % proportion of ROI sig activated
            
            recon_fold = [storage base_fold{j}(1:8) slash base_fold{j}(10:end) slash];
            flow_file = [recon_fold 'Flowc_ds.mat'];
            if exist(flow_file,'file')
                
                load(flow_file)
                flow = Fcds.signed;
                flow_roi = flow(cdata);
                pvaltmp = ind_p2tail(:,:,j,ii);
                tvaltmp = ind_t_stat(:,:,j,ii);
                flow_roi_tval = tvaltmp(cdata);
                flow_roi_pval = pvaltmp(cdata);
                flow_roi_pos_pts{i,j,ii} = find(flow_roi > 0);
                flow_roi_neg_pts{i,j,ii} = find(flow_roi < 0);
                
                flow_roi_pos_sig_pts{i,j,ii} = find(flow_roi > 1 & flow_roi_pval < Pthresh);
                flow_roi_neg_sig_pts{i,j,ii} = find(flow_roi < -1 & flow_roi_pval < Pthresh);
                
                flow_roi_pos_sig_T{i,j,ii} = nanmean(flow_roi_tval(flow_roi > 0 & flow_roi_pval < Pthresh));
                flow_roi_neg_sig_T{i,j,ii} = nanmean(flow_roi_tval(flow_roi < 0 & flow_roi_pval < Pthresh));
                
                flow_roi_pos_sig_prop{i,j,ii} = size(flow_roi_pos_sig_pts{i,j,ii},1)/size(flow_roi_pos_pts{i,j,ii},1);
                flow_roi_neg_sig_prop{i,j,ii} = size(flow_roi_neg_sig_pts{i,j,ii},1)/size(flow_roi_neg_pts{i,j,ii},1);
                
                veloc = Fcds.veloc;
                veloc_roi = veloc(cdata);
                veloc_roi_tval = tvaltmp(cdata);
                veloc_roi_pval = pvaltmp(cdata);
                veloc_roi_pos_sig_pts{i,j,ii} = find(veloc_roi > 2.5 & flow_roi_pval < Pthresh);
                veloc_roi_neg_sig_pts{i,j,ii} = find(veloc_roi < -2.5 & flow_roi_pval < Pthresh);
                
                
            end
        end
    end
    
    % Group fixed effect analysis
    grp_D = [D1 D2 D3];
    nTRs = size(grp_D,1);
    grp_c = [repmat([1 0 0 0],[1 size(base_fold,1)]), zeros(1,size(base_fold,1)), zeros(1,7*size(base_fold,1))]';
    grp_beta = zeros(size(grp_data,1),size(grp_data,2),size(grp_c,1));
    grp_t_stat = zeros(size(grp_data,1),size(grp_data,2));
    grp_p2tail = zeros(size(grp_data,1),size(grp_data,2));
    for i = 1:size(grp_data,1)
        for j = 1:size(grp_data,2)

            Y = squeeze(grp_data(i,j,:))';

            DF = nTRs - 2;

            beta_hat=inv(grp_D'*grp_D)*grp_D'*Y';
            grp_beta(i,j,:) = reshape(beta_hat,1,1,size(grp_c,1));
            Var_e=(Y'-grp_D*beta_hat)'*(Y'-grp_D*beta_hat)/DF;

            %Hypothesis testing; Compute the t statistic
            grp_t_stat(i,j)=grp_c'*beta_hat/sqrt(Var_e*grp_c'*inv(grp_D'*grp_D)*grp_c);

            if isnan(grp_t_stat(i,j))
                grp_p2tail(i,j) = 1;
            else
                grp_p2tail(i,j) = 1 - tdist2T(grp_t_stat(i,j),DF);
            end

        end
    end
    
    % Apply correction
    [grp_c_p2tail, c_alpha, h] = fwer_sidak(reshape(grp_p2tail,size(grp_p2tail,1)*size(grp_p2tail,2), []), 0.05);
    grp_c_p2tail = reshape(grp_c_p2tail,size(grp_p2tail,1),size(grp_p2tail,2));
    [h, crit_p, adj_ci_cvrg, grp_c_p2tail]=fdr_bh(grp_c_p2tail,Pthresh,'pdep','no');
%     c_p2tail = p2tail;
    
    % Average T value ROI
    for i = 1:size(roi,1)
        roifile = [roifolder roi(i).name];
        load(roifile)
        cdata = imresize(ROI.cdata,.5);

        % Average T-value also
        idx = find(cdata);
        grp_ave_t_val{i,ii} = mean(grp_t_stat(idx));

    end

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
    Cbarp = [jet(round(1.1*100))];
    h2 = figure(11);

    atmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fUS/20191119/4346075_N_D2/Anat.mat';
    anii = load(atmp); anii = anii.ref;
    % mtmp = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/fMRI/20191108/BEd_preCW_4346075_N_D5_1_47/16a/ave_anatomy_mask.nii';
    % mnii = load_nii(mtmp);
    % mask = mnii.img;
    mask = ones(size(anii));
    % Resize time-series param to anat sampling
    grp_c_p2tail_tmp = imresize(grp_c_p2tail,[size(anii,1) size(anii,2)],'bilinear');
    grp_c_p2tail_tmp(mask == 0) = 1.1;
    grp_t_stat_tmp = imresize(grp_t_stat,[size(anii,1) size(anii,2)],'bilinear');
    grp_t_stat_tmp(mask == 0) = 0;
    
    % rescale dynamic range of anat for vis
    anat = anii;
    anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
    anat = anat - maxt;
    grp_anat = anat; grp_anat(grp_c_p2tail_tmp < Pthresh) = grp_t_stat_tmp(grp_c_p2tail_tmp < Pthresh); grp_anat(mask == 0) = 0;

    figure(h2);
    subplot(size(stim,2),2,1+2*(ii-1))
    imagesc(grp_c_p2tail_tmp); set(gca,'xticklabel',[],'yticklabel',[]); colormap(gca,Cbarp); caxis([0 1.1]); title(['p-val 2tail - ' stim{ii} 'mW'])
    subplot(size(stim,2),2,2+2*(ii-1))
    imagesc(grp_anat); set(gca,'xticklabel',[],'yticklabel',[]); colormap(gca,Cbar); caxis(CAX); title(['t-val 2tail - ' stim{ii} 'mW'])

    roifolder = '/media/bradley/Seagate Backup Plus Drive/Data_Processed/segment_ofUS/fus_20191205/'; 
    roi = dir(roifolder);
    roi(contains({roi.name},'.fig')) = [];
    roi(~contains({roi.name},'R')) = [];
    
    % Time series analysis
    start = [30 69 108 147 186];
    for i = 1:size(roi,1)
        roifile = [roifolder roi(i).name];
        load(roifile)
        cdata = imresize(ROI.cdata,.5);
        
        % Total ROI
        pts{ii,i} = find(cdata);
        
        for j = 1:size(base_fold,1)
            datatmp = reshape(ind_data{j},[size(ind_data{j},1)*size(ind_data{j},2),size(ind_data{j},3)]);
            [b,a] = butter(4,[0.0005 0.1]/((1/1.5)/2));
            
            ts{ii,i}(j,:) = sum(datatmp(pts{ii,i},:));
            meanbase = mean(ts{ii,i}(j,1:20));
            tsnorm{ii,i}(j,:) = (ts{ii,i}(j,:) - meanbase)/meanbase*100;
            tsnorm{ii,i}(j,:) = filtfilt(b,a,tsnorm{ii,i}(j,:));
            
            ts_normpeaks = zeros(size(start,2)-1,size(-10:28,2));
            for k = 1:size(start,2)
                ts_normpeaks(k,:) = tsnorm{ii,i}(j,start(k)-10:start(k)+28);
            end
            tmp = ts_normpeaks'-repmat(mean(ts_normpeaks',1),[39 1]);
            tmp = tmp - repmat(mean(tmp(1:10,:),1),[39,1]);
            tsnormave{ii,i}(j,:) = mean(tmp',1);
            tsnormavepeak{ii,i}(j) = max(mean(tmp'));
            tsnormaveAUC{ii,i}(j) = sum(mean(tmp'));
            
        end
        
        % Vascular contributions to ROI
        pts_pos = flow_roi_pos_pts{i,j,ii};
        pts_neg = flow_roi_neg_pts{i,j,ii};
        
        pts_pos_act = flow_roi_pos_sig_pts{i,j,ii};
        pts_neg_act = flow_roi_neg_sig_pts{i,j,ii};
        
        for j = 1:size(base_fold,1)
            datatmp = reshape(ind_data{j},[size(ind_data{j},1)*size(ind_data{j},2),size(ind_data{j},3)]);
            [b,a] = butter(4,[0.01 .4]);
            
            ts_tot = datatmp(pts{ii,i},:);
            ts_pos = sum(ts_tot(flow_roi_pos_pts{i,j,ii},:),1);
            ts_neg = sum(ts_tot(flow_roi_neg_pts{i,j,ii},:),1);
            ts_sig_pos = sum(ts_tot(flow_roi_pos_sig_pts{i,j,ii},:),1);
            ts_sig_neg = sum(ts_tot(flow_roi_neg_sig_pts{i,j,ii},:),1);
            ts_sig_fast = sum(ts_tot(veloc_roi_pos_sig_pts{i,j,ii},:),1);
            ts_sig_slow = sum(ts_tot(veloc_roi_neg_sig_pts{i,j,ii},:),1);
            
            meanbase_pos = mean(ts_pos(1:20));
            meanbase_neg = mean(ts_neg(1:20));
            meanbase_sig_pos = mean(ts_sig_pos(1:20));
            meanbase_sig_neg = mean(ts_sig_neg(1:20));
            meanbase_sig_fast = mean(ts_sig_fast(1:20));
            meanbase_sig_slow = mean(ts_sig_slow(1:20));
            
            NumA = size(flow_roi_pos_sig_pts{i,j,ii},1);
            NumN = size(flow_roi_neg_sig_pts{i,j,ii},1);
            Aprop{ii,i}(j) = NumA;% / (NumA + NumN) *100;
            Nprop{ii,i}(j) = NumN;% / (NumA + NumN) *100;
            
            tsnorm_pos{ii,i}(j,:) = (ts_pos - meanbase_pos)/meanbase_pos*100;
            tsnorm_neg{ii,i}(j,:) = (ts_neg - meanbase_neg)/meanbase_neg*100;
            tsnorm_sig_pos{ii,i}(j,:) = (ts_sig_pos - meanbase_sig_pos)/meanbase_sig_pos*100;
            tsnorm_sig_neg{ii,i}(j,:) = (ts_sig_neg - meanbase_sig_neg)/meanbase_sig_neg*100;
            tsnorm_sig_fast{ii,i}(j,:) = (ts_sig_fast - meanbase_sig_fast)/meanbase_sig_fast*100;
            tsnorm_sig_slow{ii,i}(j,:) = (ts_sig_slow - meanbase_sig_slow)/meanbase_sig_slow*100;
            
            tsnorm_pos{ii,i}(j,:) = filtfilt(b,a,tsnorm_pos{ii,i}(j,:));
            tsnorm_neg{ii,i}(j,:) = filtfilt(b,a,tsnorm_neg{ii,i}(j,:));
            tsnorm_sig_pos{ii,i}(j,:) = filtfilt(b,a,tsnorm_sig_pos{ii,i}(j,:));
            tsnorm_sig_neg{ii,i}(j,:) = filtfilt(b,a,tsnorm_sig_neg{ii,i}(j,:));
            tsnorm_sig_fast{ii,i}(j,:) = filtfilt(b,a,tsnorm_sig_fast{ii,i}(j,:));
            tsnorm_sig_slow{ii,i}(j,:) = filtfilt(b,a,tsnorm_sig_slow{ii,i}(j,:));
            
            ts_normpeaks_pos = zeros(size(start,2),size(-10:28,2));
            ts_normpeaks_neg = zeros(size(start,2),size(-10:28,2));
            ts_normpeaks_sig_pos = zeros(size(start,2),size(-10:28,2));
            ts_normpeaks_sig_neg = zeros(size(start,2),size(-10:28,2));
            ts_normpeaks_sig_fast = zeros(size(start,2),size(-10:28,2));
            ts_normpeaks_sig_slow = zeros(size(start,2),size(-10:28,2));
            for k = 1:size(start,2)
                ts_normpeaks_pos(k,:) = tsnorm_pos{ii,i}(j,start(k)-10:start(k)+28);
                ts_normpeaks_neg(k,:) = tsnorm_neg{ii,i}(j,start(k)-10:start(k)+28);
                ts_normpeaks_sig_pos(k,:) = tsnorm_sig_pos{ii,i}(j,start(k)-10:start(k)+28);
                ts_normpeaks_sig_neg(k,:) = tsnorm_sig_neg{ii,i}(j,start(k)-10:start(k)+28);
                ts_normpeaks_sig_fast(k,:) = tsnorm_sig_fast{ii,i}(j,start(k)-10:start(k)+28);
                ts_normpeaks_sig_slow(k,:) = tsnorm_sig_slow{ii,i}(j,start(k)-10:start(k)+28);
            end
            tmp1 = ts_normpeaks_sig_pos'-repmat(mean(ts_normpeaks_sig_pos',1),[39 1]);
            tmp1 = tmp1 - repmat(mean(tmp1(1:10,:),1),[39,1]);
            tmp2 = ts_normpeaks_sig_neg'-repmat(mean(ts_normpeaks_sig_neg',1),[39 1]);
            tmp2 = tmp2 - repmat(mean(tmp2(1:10,:),1),[39,1]);
            tmp3 = ts_normpeaks_sig_fast'-repmat(mean(ts_normpeaks_sig_fast',1),[39 1]);
            tmp3 = tmp3 - repmat(mean(tmp3(1:10,:),1),[39,1]);
            tmp4 = ts_normpeaks_sig_slow'-repmat(mean(ts_normpeaks_sig_slow',1),[39 1]);
            tmp4 = tmp4 - repmat(mean(tmp4(1:10,:),1),[39,1]);
            tsnormave_pos{ii,i}(j,:) = mean(ts_normpeaks_pos,1);
            tsnormave_neg{ii,i}(j,:) = mean(ts_normpeaks_neg,1);
            tsnormave_sig_pos{ii,i}(j,:) = mean(tmp1',1);
            tsnormave_sig_neg{ii,i}(j,:) = mean(tmp2',1);
            tsnormave_sig_fast{ii,i}(j,:) = mean(tmp3',1);
            tsnormave_sig_slow{ii,i}(j,:) = mean(tmp4',1);
            
            tsnormavepeak_pos(i,j,ii) = mean(max(tsnormave_pos{ii,i}(j,:)'));
            tsnormavepeak_neg(i,j,ii) = mean(max(tsnormave_neg{ii,i}(j,:)'));
            tsnormavepeak_sig_pos(i,j,ii) = max(mean(tmp1'));
            tsnormavepeak_sig_neg(i,j,ii) = max(mean(tmp2'));
            tsnormavepeak_sig_fast(i,j,ii) = max(mean(tmp3'));
            tsnormavepeak_sig_slow(i,j,ii) = max(mean(tmp4'));
            tsnormaveAUC_sig_pos(i,j,ii) = sum(mean(tmp1,2));
            tsnormaveAUC_sig_neg(i,j,ii) = sum(mean(tmp2,2));
            tsnormaveAUC_sig_fast(i,j,ii) = sum(mean(tmp3,2));
            tsnormaveAUC_sig_slow(i,j,ii) = sum(mean(tmp4,2));
            
        end
        
        [H,PP1(i,ii)] = ttest(tsnormavepeak_pos(i,:,ii),tsnormavepeak_neg(i,:,ii));
%         [H,PP2(i,ii)] = ttest(tsnormavepeak_sig_pos(i,:,ii),tsnormavepeak_sig_neg(i,:,ii));
    end
    
    % ROI Analysis
    maxp = 125;
    anattmp = -2*maxp*ones(size(imresize(anat,2)));
    for i = 1:size(roi,1)
        roifile = [roifolder roi(i).name];
        load(roifile)
        cdata = imresize(ROI.cdata,2);
        pts2 {i}= find(cdata);
        anattmp(pts2{i}) = mean(tsnormavepeak{ii,i},2);
    end
%     maxp = max(max(anattmp));
    roi_anat = imresize(anat + maxt,2);
    roi_anat = roi_anat - maxp;
    roi_anat(anattmp ~= -2*maxp) = anattmp(anattmp ~= -2*maxp);
    
    Cbar = [gray(abs(plotmax - plotmin)*100); jet(2*maxp*100)]; CAX = [plotmin - maxp maxp];
    figure; subplot(1,2,1); imagesc(roi_anat); caxis(CAX);
    subplot(1,2,2); imagesc(anattmp);colormap(Cbar); caxis(CAX);

end

% Save a few variables: time series peak, AUC
save([storage 'fUSI_time_series_peak_' descr '.mat'],'tsnormavepeak','fusroi');
save([storage 'fUSI_time_series_AUC_' descr '.mat'],'tsnormaveAUC','fusroi');

save([storage 'fUSI_time_series_peak_ven_' descr '.mat'],'tsnormavepeak_sig_neg','fusroi');
save([storage 'fUSI_time_series_peak_art_' descr '.mat'],'tsnormavepeak_sig_pos','fusroi');


% Average time series across ROIs and intensities
figure(500); clf
for j = 1:20
    subplot(5,4,j); hold on

    tmp3 = tsnormave{3,j};
    stdshade(tmp3,0.25,[32 104 52]/255)
    tmp2 = tsnormave{2,j};
    stdshade(tmp2,0.25,[90 184 70]/255)
    tmp1 = tsnormave{1,j};
    stdshade(tmp1,0.25,[175 216 166]/255)
    
    if strcmp(descr,'ctrl')
        set(gca,'ylim',[-2.5 10])
    else
        set(gca,'ylim',[-50 150])
    end
    y1=get(gca,'ylim');
    plot([11 11],y1,'k'); plot([19 19],y1,'k')
    title(fusroi{j})
end
suptitle('ROI Time Series: Intensity');

% VASCULAR INFO
figure(501); clf; figure(502); clf; figure(503); clf; figure(504); clf;
figure(505); clf; figure(506); clf; figure(507); clf; figure(508); clf
clear StatAUC StatPk StatCount
StatAUC(1,1:5) = {'subjid','vessel','intensity','pairInt','pairVes'};
id = repmat(1:size(base_fold,1),3,1); id = repmat(id(:),2,1);
StatAUC(2:1+size(base_fold,1)*6,1) = num2cell(id);
StatAUC(2:1+size(base_fold,1)*3,2) = {'Art'};
StatAUC(2+size(base_fold,1)*3:1+size(base_fold,1)*6,2) = {'Vein'};
StatAUC(2:1+size(base_fold,1)*6,3) = num2cell(repmat([1;2;3],2*size(base_fold,1),1));
StatAUC(2:1+size(base_fold,1)*6,4) = num2cell(repmat([1;2;3],2*size(base_fold,1),1));
StatAUC(2:1+size(base_fold,1)*3,5) = num2cell(repmat([1;2;3],size(base_fold,1),1));
StatAUC(2+size(base_fold,1)*3:1+size(base_fold,1)*6,5) = num2cell(repmat([4;5;6],size(base_fold,1),1));
StatPk = StatAUC; StatCount = StatAUC;

for j = 1:20
    % Plot vascular TS - 0.1 mW
    figure(501); subplot(5,4,j); hold on
    tmp1 = tsnormave_sig_pos{1,j}; stdshade(tmp1,0.25,'b')
    tmp2 = tsnormave_sig_neg{1,j}; stdshade(tmp2,0.25,'r')
    if max([max(mean(tmp1,1)) max(mean(tmp2,1))]) > 50; M = 150; else; M = 75; end
    set(gca,'ylim',[-25 M]); y1=get(gca,'ylim');
    plot([11 11],y1,'k'); plot([19 19],y1,'k'); title(fusroi{j})
    
    % Plot vascular TS - 0.5 mW
    figure(502); subplot(5,4,j); hold on
    tmp1 = tsnormave_sig_pos{2,j}; stdshade(tmp1,0.25,'b')
    tmp2 = tsnormave_sig_neg{2,j}; stdshade(tmp2,0.25,'r')
    if max([max(mean(tmp1,1)) max(mean(tmp2,1))]) > 50; M = 150; else; M = 75; end
    set(gca,'ylim',[-25 M]); y1=get(gca,'ylim');
    plot([11 11],y1,'k'); plot([19 19],y1,'k'); title(fusroi{j})
    
    % Plot vascular TS - 1.0 mW
    figure(503); subplot(5,4,j); hold on
    tmp1 = tsnormave_sig_pos{3,j}; stdshade(tmp1,0.25,'b')
    tmp2 = tsnormave_sig_neg{3,j}; stdshade(tmp2,0.25,'r')
    if max([max(mean(tmp1,1)) max(mean(tmp2,1))]) > 50; M = 150; else; M = 75; end
    set(gca,'ylim',[-25 M]); y1=get(gca,'ylim');
    plot([11 11],y1,'k'); plot([19 19],y1,'k'); title(fusroi{j})
    
    % Plot Arteriole TS across all intensities
    figure(504); subplot(5,4,j); hold on
    tmp3 = tsnormave_sig_pos{3,j}; stdshade(tmp3,0.25,'k')
    tmp2 = tsnormave_sig_pos{2,j}; stdshade(tmp2,0.25,'m')
    tmp1 = tsnormave_sig_pos{1,j}; stdshade(tmp1,0.25,'c')
    
    if max([max(mean(tmp1,1)) max(mean(tmp2,1)) max(mean(tmp3,1))]) > 50; M = 150; else; M = 75; end
    set(gca,'ylim',[-50 M],'xlim',[11 30]); y1=get(gca,'ylim');
    plot([11 11],y1,'k'); plot([19 19],y1,'k'); title(fusroi{j})
    
    % Plot Venule TS across all intensities
    figure(505); subplot(5,4,j); hold on
    tmp3 = tsnormave_sig_neg{3,j}; stdshade(tmp3,0.25,'k')
    tmp2 = tsnormave_sig_neg{2,j}; stdshade(tmp2,0.25,'m')
    tmp1 = tsnormave_sig_neg{1,j}; stdshade(tmp1,0.25,'c')
    if max([max(mean(tmp1,1)) max(mean(tmp2,1)) max(mean(tmp3,1))]) > 50; M = 150; else; M = 75; end
    set(gca,'ylim',[-50 M],'xlim',[11 30]); y1=get(gca,'ylim');
    plot([11 11],y1,'k'); plot([19 19],y1,'k'); title(fusroi{j})
    
    % Plot TS Area Under Curve for vascular components
    figure(506); subplot(5,4,j); hold on
    AUC_pos_mean = squeeze(nanmean(tsnormaveAUC_sig_pos(j,:,:),2));
    AUC_pos_std = squeeze(nanstd(tsnormaveAUC_sig_pos(j,:,:),0))/sqrt(size(base_fold,1));
    AUC_neg_mean = squeeze(nanmean(tsnormaveAUC_sig_neg(j,:,:),2));
    AUC_neg_std = squeeze(nanstd(tsnormaveAUC_sig_neg(j,:,:),0))/sqrt(size(base_fold,1));
    Ave = [AUC_pos_mean AUC_neg_mean]; Sem = [AUC_pos_std AUC_neg_std];
    ctrs = 1:2; hBar{1} = bar(ctrs, Ave'); drawnow; clear ctr ydt
    for k1 = 1:size(Ave',2)
        ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
        ydt(k1,:) = hBar{1}(k1).YData;
    end
    errorbar(ctr, ydt, zeros(2,3)', Sem, zeros(2,3)', zeros(2,3)', '.k',...
        'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
    hBar{2} = bar(ctrs, Ave'); delete(hBar{1})
    if max(Ave(:)) < 750; M = 1000; else M = 2000; end
    set(gca,'ylim',[0 M],'xtick',[1 2],'xticklabel',{'art' 'ven'});
    test_pos = []; test_neg = [];
    for k = 1:3
        tmp1 = tsnormaveAUC_sig_pos(j,:,k); tmp11 = tmp1(~isnan(tmp1));
        tmp2 = tsnormaveAUC_sig_neg(j,:,k); tmp22 = tmp2(~isnan(tmp2));
        if size(tmp11,2) == size(tmp22,2) && size(tmp11,2)>1
            [H,P] = ttest(tmp11,tmp22);
%             [P,H] = signrank(tmp1,tmp2);
            text(ctr(k)-0.15,max(Ave(:))+max(Sem(:))+100,num2str(P,'%1.2f'))
        end
        test_pos = [test_pos tmp1']; test_neg = [test_neg tmp2'];
    end
    
    % Formate values for R txt file
    StatAUC{1,5+j} = fusroi{j};
    test_posR = []; test_negR = [];
    for k = 1:size(base_fold,1)
        test_posR = [test_posR; test_pos(k,:)'];
        test_negR = [test_negR; test_neg(k,:)'];
    end
    StatAUC(2:1+size(base_fold,1)*3,5+j) = num2cell(test_posR);
    StatAUC(2+size(base_fold,1)*3:1+size(base_fold,1)*6,5+j) = num2cell(test_negR);
%     set(gca,'ylim',[0 max(Ave(:))+max(Sem(:))+200],'xtick',[]);
    title(fusroi{j});hold off
    
    % Plot TS peak for vascular components
    figure(507); subplot(5,4,j); cla; hold on
    Pk_pos_mean = squeeze(nanmean(tsnormavepeak_sig_pos(j,:,:),2));
    Pk_pos_std = squeeze(nanstd(tsnormavepeak_sig_pos(j,:,:),0))/sqrt(size(base_fold,1));
    Pk_neg_mean = squeeze(nanmean(tsnormavepeak_sig_neg(j,:,:),2));
    Pk_neg_std = squeeze(nanstd(tsnormavepeak_sig_neg(j,:,:),0))/sqrt(size(base_fold,1));
    Ave = [Pk_pos_mean Pk_neg_mean]; Sem = [Pk_pos_std Pk_neg_std];
    ctrs = 1:2; hBar{1} = bar(ctrs, Ave'); drawnow; clear ctr ydt
    for k1 = 1:size(Ave',2)
        ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
        ydt(k1,:) = hBar{1}(k1).YData;
    end
    errorbar(ctr, ydt, zeros(2,3)', Sem, zeros(2,3)', zeros(2,3)', '.k',...
        'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
    hBar{2} = bar(ctrs, Ave'); delete(hBar{1})
    test_pos = []; test_neg = [];
    for k = 1:3
        tmp1 = tsnormavepeak_sig_pos(j,:,k); tmp11 = tmp1(~isnan(tmp1));
        tmp2 = tsnormavepeak_sig_neg(j,:,k); tmp22 = tmp2(~isnan(tmp2));
        if size(tmp11,2) == size(tmp22,2) && size(tmp11,2)>1
            [H,P] = ttest(tmp11,tmp22);
%             [P,H] = signrank(tmp1,tmp2);
            text(ctr(k)-0.15,max(Ave(:))+max(Sem(:))+100,num2str(P,'%1.2f'))
        end
        test_pos = [test_pos tmp1']; test_neg = [test_neg tmp2'];
    end
    
    % Formate values for R txt file
    StatPk{1,5+j} = fusroi{j};
    test_posR = []; test_negR = [];
    for k = 1:size(base_fold,1)
        test_posR = [test_posR; test_pos(k,:)'];
        test_negR = [test_negR; test_neg(k,:)'];
    end
    StatPk(2:1+size(base_fold,1)*3,5+j) = num2cell(test_posR);
    StatPk(2+size(base_fold,1)*3:1+size(base_fold,1)*6,5+j) = num2cell(test_negR);
%     set(gca,'ylim',[0 max(Ave(:))+max(Sem(:))+200],'xtick',[]); title(fusroi{j})
    hold off
    
    % Plot TS Area Under Curve for flow speeds
    figure(508); subplot(5,4,j); hold on
    AUC_fast_mean = squeeze(nanmean(tsnormaveAUC_sig_fast(j,:,:),2));
    AUC_fast_std = squeeze(nanstd(tsnormaveAUC_sig_fast(j,:,:),0))/sqrt(size(base_fold,1));
    AUC_slow_mean = squeeze(nanmean(tsnormaveAUC_sig_slow(j,:,:),2));
    AUC_slow_std = squeeze(nanstd(tsnormaveAUC_sig_slow(j,:,:),0))/sqrt(size(base_fold,1));
    Ave = [AUC_fast_mean AUC_slow_mean]; Sem = [AUC_fast_std AUC_slow_std];
    ctrs = 1:2; hBar{1} = bar(ctrs, Ave'); drawnow; clear ctr ydt
    for k1 = 1:size(Ave',2)
        ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
        ydt(k1,:) = hBar{1}(k1).YData;
    end
    errorbar(ctr, ydt, zeros(2,3)', Sem, zeros(2,3)', zeros(2,3)', '.k',...
        'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
    hBar{2} = bar(ctrs, Ave'); delete(hBar{1})
    if max(Ave(:)) < 750; M = 1000; else M = 5000; end
    set(gca,'ylim',[0 M],'xtick',[1 2],'xticklabel',{'fast' 'slow'}); title(fusroi{j})
    test_fast = []; test_slow = [];
    for k = 1:3
        tmp1 = tsnormaveAUC_sig_fast(j,:,k); tmp11 = tmp1(~isnan(tmp1));
        tmp2 = tsnormaveAUC_sig_slow(j,:,k); tmp22 = tmp2(~isnan(tmp2));
        if size(tmp11,2) == size(tmp22,2) && size(tmp11,2)>1
            [H,P] = ttest(tmp11,tmp22);
%             [P,H] = signrank(tmp1,tmp2);
            text(ctr(k)-0.15,max(Ave(:))+max(Sem(:))+100,num2str(P,'%1.2f'))
        end
        test_fast = [test_fast tmp1']; test_slow = [test_slow tmp2'];
    end

end
figure(501); suptitle('0.1mW V vs A');
figure(502); suptitle('0.5mW V vs A');
figure(503); suptitle('1.0mW V vs A');
figure(504); suptitle('Arterioles TS');
figure(505); suptitle('Venules TS');
figure(506); suptitle('Vascular AUC');
figure(507); suptitle('Vascular Peak');
figure(508); suptitle('Vascular Speed');

% Write stat file for R - AUC vascular components
StatAUC(:,[11,25]) = []; % Remove S2 regions!
fid = fopen(['/home/bradley/Dropbox/Vascular_AUC_' descr '.txt'],'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatAUC{1,:});
for K = 2:size(StatAUC,1)
    fprintf(fid, '%.0f %s %.2f %.0f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatAUC{K,:});
end
fclose(fid)

% Write stat file for R - Peak Response vascular components
StatPk(:,[11,25]) = [];
fid = fopen(['/home/bradley/Dropbox/Vascular_Pk_' descr '.txt'],'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatPk{1,:});
for K = 2:size(StatPk,1)
    fprintf(fid, '%.0f %s %.2f %.0f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatPk{K,:});
end
fclose(fid)

% Pixel count of activated ROI Arteriole vs Venule
figure(509); clf
for i = 1:20
    for j = 1:3
        meanAprop(j,i) = mean(Aprop{j,i});
        meanNprop(j,i) = mean(Nprop{j,i});
        Ta = Aprop{j,i}; Tn = Nprop{j,i};
        N = [Ta' Tn']; Nidx = sum(isnan(N),2);
        Ta(find(Nidx > 0)) = []; Tn(find(Nidx > 0)) = [];
        if ~isempty(Ta) && ~isempty(Tn) && size(Tn,2) > 1
            [H,P(j,i)] = ttest(Ta,Tn);
            [P(j,i),H] = signrank(Ta,Tn);
        else
            P(j,i) = 1;
        end
    end
    
    subplot(5,4,i); cla; hold on
    Ave = [nanmean(vertcat(Aprop{:,i}),2) nanmean(vertcat(Nprop{:,i}),2)];
    Sem = [nanstd(vertcat(Aprop{:,i}),0,2)/sqrt(7) nanstd(vertcat(Nprop{:,i}),0,2)/sqrt(7)];
    ctrs = 1:2; hBar{1} = bar(ctrs, Ave'); drawnow; clear ctr ydt
    for k1 = 1:size(Ave',2)
        ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
        ydt(k1,:) = hBar{1}(k1).YData;
    end
    errorbar(ctr, ydt, zeros(2,3)', Sem, zeros(2,3)', zeros(2,3)', '.k',...
        'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
    hBar{2} = bar(ctrs, Ave'); delete(hBar{1})
    if max(Ave(:)) > 100; M = 350; else; M = 150; end
    set(gca,'ylim',[0 M],'xtick',[1 2],'xticklabel',{'art' 'ven'}); title(fusroi{i})
    
    
    StatCount{1,5+i} = fusroi{i};
    StatCount(2:1+size(base_fold,1)*3,5+i) = num2cell(horzcat(Aprop{:,i})');
    StatCount(2+size(base_fold,1)*3:1+size(base_fold,1)*6,5+i) = num2cell(horzcat(Nprop{:,i})');
    
end
suptitle('Vascular Component Count');

% Write stat file for R - Vascular count
StatCount(:,[11,25]) = [];
fid = fopen(['/home/bradley/Dropbox/Vascular_Count_' descr '.txt'],'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatCount{1,:});
for K = 2:size(StatCount,1)
    fprintf(fid, '%.0f %s %.2f %.0f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatCount{K,:});
end
fclose(fid)



% Proportion ROI activated
Prop = cell2mat(ind_p2_roi_prop);
PropM = squeeze(mean(Prop,2));
PropSem = squeeze(std(Prop,0,2)/sqrt(size(Prop,2)));

figure;
group = cell(size(base_fold,1),1); group(:) = {'1'};
for i = 1:size(roi,1)
    for j = 1:size(base_fold,1)
        Prop2(j,:,i) = squeeze(Prop(i,j,:))';
    end
    
    subplot(5,4,i); hold on
    parallelcoords(Prop2(:,:,i),'Group',group,'Labels',stim,'color','k');
    errorbar(1:3,mean(Prop2(:,:,i),1),std(Prop2(:,:,i),1)/sqrt(size(Prop2(:,:,i),1)),...
        '-o','markersize',3,'linewidth',1.5,'color','r')
    title(fusroi{i})

    legend off; set(gca,'ylabel',[],'ylim',[0 1.25],'xlim',[0.75 3.25])
    
    [H,P12] = ttest(Prop2(:,1,i),Prop2(:,2,i)); text(1.35,1,num2str(P12,'%1.2f'))
    [H,P23] = ttest(Prop2(:,2,i),Prop2(:,3,i)); text(2.35,1,num2str(P23,'%1.2f'))
    [H,P13] = ttest(Prop2(:,1,i),Prop2(:,3,i)); text(1.85,1.25,num2str(P13,'%1.2f'))
end
suptitle('Proportion ROI Activated')

ctrs = 1:20;
data = PropM;
figure; hBar{1} = bar(ctrs, data); hold on
clear ctr ydt
for k1 = 1:size(PropM,2)
    ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
    ydt(k1,:) = hBar{1}(k1).YData;
end
errorbar(ctr, ydt, zeros(3,20), PropSem', zeros(3,20), zeros(3,20), '.k',...
    'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
hBar{2} = bar(ctrs, data); delete(hBar{1})
set(gca,'ylim',[0 1],'xtick',[]); hold off
suptitle('Proportion ROI Activated')

% Ave t-stat
Tval = cell2mat(ind_t_stat_ave);
TvalM = squeeze(nanmean(Tval,2));
TvalSem = squeeze(nanstd(Tval,0,2)/sqrt(size(Tval,2)));
ctrs = 1:20;
data = TvalM;
figure; hBar{1} = bar(ctrs, data); hold on
clear ctr ydt
for k1 = 1:size(TvalM,2)
    ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
    ydt(k1,:) = hBar{1}(k1).YData;
end
errorbar(ctr, ydt, zeros(3,20), TvalSem', zeros(3,20), zeros(3,20), '.k',...
    'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
hBar{2} = bar(ctrs, data); delete(hBar{1})
set(gca,'ylim',[0 11],'xtick',1:20,'xticklabel',{roi.name})
hold off; suptitle('Average T-value')

clear StatROIt
StatROIt(1,1:3) = {'subjid','intensity','pairInt'};
id = repmat(1:size(base_fold,1),3,1); StatROIt(2:1+size(base_fold,1)*3,1) = num2cell(id(:));
StatROIt(2:1+size(base_fold,1)*3,2) = num2cell(repmat([1;2;3],size(base_fold,1),1));
StatROIt(2:1+size(base_fold,1)*3,3) = num2cell(repmat([1;2;3],size(base_fold,1),1));
for  i = 1:20
        tmp1 = cell2mat(squeeze(ind_t_stat_ave(i,:,:)))';
        StatROIt(1,3+i) = fusroi(i);
        StatROIt(2:1+size(base_fold,1)*3,3+i) = num2cell(tmp1(:));
end

% Write stat file for R - ROI average T values
StatROIt(:,[9,23]) = [];
fid = fopen(['/home/bradley/Dropbox/ROI_Tval_fUSI_' descr '.txt'],'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatROIt{1,:});
for K = 2:size(StatROIt,1)
    fprintf(fid, '%.0f %.2f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatROIt{K,:});
end
fclose(fid)

%%

%%%%%%%% REGRESSION
BB = zeros(size(ind_t_stat,1),size(ind_t_stat,2));
RR = zeros(size(ind_t_stat,1),size(ind_t_stat,2));
PP = zeros(size(ind_t_stat,1),size(ind_t_stat,2));
for i = 1:size(ind_t_stat,1)
    for j = 1:size(ind_t_stat,2)
        
        tmp1 = squeeze(ind_t_stat(i,j,:,1));
        tmp2 = squeeze(ind_t_stat(i,j,:,2));
        tmp3 = squeeze(ind_t_stat(i,j,:,3));
        
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
        
        if i == 10 && j == 46
            
            figure; scatter(T2,T1); hold on
            set(gca,'xlim',[-.5 1.5],'ylim',[-5 35]);
            plot(reshape(T2,[7,3])',reshape(T1,[7,3])')
            scatter([0.1 0.5 1.0], mean(reshape(T1,[7,3]),1),50,'k','filled')
            plot([0.1 0.5 1.0],mean(reshape(T1,[7,3]),1),'k','linewidth',2)
            
        end
        
    end
end


% rescale dynamic range of anat for vis
maxR = max(RR(:)); maxB = ceil(max(abs(BB(:)))); plotmin = -100; plotmax = 0;
CbarR = [gray(abs(plotmax - plotmin)*100); fireice(2*maxR*100)]; CAXR = [plotmin - maxR maxR];
CbarB = [gray(abs(plotmax - plotmin)*100); fireice(2*maxB*100)]; CAXB = [plotmin - maxB maxB];

BB2 = imresize(BB,size(anii));
RR2 = imresize(RR,size(anii));
PP2 = imresize(PP,size(anii));
anat = anii;
anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin; 
anatRR = anat - maxR; anatBB = anat - maxB;
anatBB(PP2 < .05) = BB2(PP2 < .05); anatBB(mask(:,:) == 0) = CAXB(1);
anatRR(PP2 < .05) = RR2(PP2 < .05); anatRR(mask(:,:) == 0) = CAXR(1);

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



