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
storage = [Dir 'Data_Processed' slash 'ephys' slash]; if ~exist(storage,'dir'); mkdir(storage); end
addpath(genpath(storage));
raw = [Dir 'Ephys'];

% Total Thy1
base_fold = {'20191217_4346075_N';
    '20191217_4364123_R1';
    '20191217_4364122_N';
    '20191217_4364143_R1';
    };
type = 'Thy1';

% Control
% base_fold = {'20200302_4419409_N';
%     '20200302_4419410_N';
%     };
% type = 'ctrl';

%%
% Read in raw data and resave condensed versions
final_GetPlxData_BE(storage,base_fold,slash,type)
% Break data into trials for LFP and spiking
final_RearrangeData_BE(storage,base_fold,slash);
% Check clippings
final_CheckLFPtraces_BE(storage,base_fold,slash);
% Quality control 1
final_RunQC1_BE(storage,base_fold,slash);
% Quality control 2
final_RunQC2_BE(storage,base_fold,slash);
% Drop bad channels and trials detected during QC2
final_CleanData_BE(storage,base_fold,slash,type);
% Check spiking activity (if clipping affects behavior)
final_CheckSpikesDuringClips_BE(storage,base_fold,slash);

% LFP
% plot and save average LFP for each channel and region (individual animal)
final_PlotLFPspectrograms_Averaged_BE(storage,base_fold,slash);
% Calculate LFP amplitude (10 Hz)
final_GetLFPamplitude_Averaged_BE(storage,base_fold,slash)
% Group LFP Trace and Spectrogram Plots
final_PlotLFP_Info_BE(storage,base_fold,slash,type);

% SPIKES
% Plot and save individual spike waveforms
final_CheckWaveforms_BE(storage,base_fold,slash)
% Significant changes to spiking
final_TestSpikes_BE(storage,base_fold,slash);
% Plot raster plots for different neurons for each animal/stim
final_PlotRaster_BE(storage,base_fold,slash);
% Plot stacked bar graphs for significant changes during stim (group)
final_PlotStackedBar_BE(storage,base_fold,slash)
% Plot firing rate per time bin
final_PlotFiringRatesBin_BE(storage,base_fold,slash)
% Plot firing rate per stimulus period
final_PlotFiringRatesPeriod_BE(storage,base_fold,slash)

%%
% Some neurovascular coupling analysis

stim = {'0_1' '0_5' '1_0'};
for i_stim = 1:size(stim,2)
    
    storage = [Dir 'Data_Processed' slash 'fUS' slash];
    stim_storage = [storage stim{i_stim} slash 'ephys' slash];
    stim_file = [stim_storage stim{i_stim} '_time_series_info.mat'];
    load(stim_file)
    fusi(:,:,i_stim) = cell2mat(tsnormaveZ_roi);
    fusi_noise(:,:,i_stim) = cell2mat(tsnormaveZ_noise);
    
    storage = [Dir 'Data_Processed' slash 'fMRI' slash];
    stim_storage = [storage stim{i_stim} slash 'ephys' slash];
    stim_file = [stim_storage stim{i_stim} '_time_series_info.mat'];
    load(stim_file)
    fmri(:,:,i_stim) = cell2mat(tsnormaveZ_roi);
    fmri_noise(:,:,i_stim) = cell2mat(tsnormaveZ_noise);
    
    % controls
    storage = [Dir 'Data_Processed' slash 'fUS' slash];
    stim_storage = [storage stim{i_stim} slash 'ctrl' slash];
    stim_file = [stim_storage stim{i_stim} '_time_series_info.mat'];
    load(stim_file)
    fusi_ctrl(:,:,i_stim) = cell2mat(tsnormaveZ_roi);
    
    storage = [Dir 'Data_Processed' slash 'fMRI' slash];
    stim_storage = [storage stim{i_stim} slash 'ctrl' slash];
    stim_file = [stim_storage stim{i_stim} '_time_series_info.mat'];
    load(stim_file)
    fmri_ctrl(:,:,i_stim) = cell2mat(tsnormaveZ_roi);

end

%%
% NOISE/CTRL Analysis
ridx = [1 5 14 15];
T = {'lcpu' 'lm1' 'rm1' 'rcpu'};
for i_region = 1:size(ridx,2)
    for i_stim = 1:3
        FMRI{i_region,i_stim} = fmri(ridx(i_region),:,i_stim);
        FMRI_noise{i_region,i_stim} = fmri_noise(:,:,i_stim);
        FMRI_ctrl{i_region,i_stim} = fmri_ctrl(ridx(i_region),:,i_stim);
        
        
        FUSI{i_region,i_stim} = fusi(ridx(i_region),:,i_stim);
        FUSI_noise{i_region,i_stim} = fusi_noise(:,:,i_stim);
        FUSI_ctrl{i_region,i_stim} = fusi_ctrl(ridx(i_region),:,i_stim);
        
        figure; subplot(1,2,1); bar([mean(FMRI{i_region,i_stim}) mean(FMRI_noise{i_region,i_stim}) mean(FMRI_ctrl{i_region,i_stim})]);
        hold on; errorbar([mean(FMRI{i_region,i_stim}) mean(FMRI_noise{i_region,i_stim}) mean(FMRI_ctrl{i_region,i_stim})],...
            [std(FMRI{i_region,i_stim}) std(FMRI_noise{i_region,i_stim}) std(FMRI_ctrl{i_region,i_stim})]/2)
        subplot(1,2,2); bar([mean(FUSI{i_region,i_stim}) mean(FUSI_noise{i_region,i_stim}) mean(FUSI_ctrl{i_region,i_stim})]);
        hold on; errorbar([mean(FUSI{i_region,i_stim}) mean(FUSI_noise{i_region,i_stim}) mean(FUSI_ctrl{i_region,i_stim})],...
            [std(FUSI{i_region,i_stim}) std(FUSI_noise{i_region,i_stim}) std(FUSI_ctrl{i_region,i_stim})]/2)
        
    end
end

figure;
Zscore = {'subjid' 'mod' 'intensity' ' pairInt' 'pairMod'};
Zscore(2:25,1) = num2cell(repmat([1:4],1,6)');
Zscore(2:13,2) = repmat({'FMRI'},12,1); Zscore(14:25,2) = repmat({'FUSI'},12,1);
Zscore(2:25,3) = num2cell(repmat([1 1 1 1 2 2 2 2 3 3 3 3],1,2)');
Zscore(2:25,4) = num2cell(repmat([1 1 1 1 2 2 2 2 3 3 3 3],1,2)');
Zscore(2:25,5) = num2cell([1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]);

for i_region = 1:4
    figure(25); subplot(2,2,i_region); hold on
    errorbar([0.1 .5 1], [mean(FMRI{i_region,1}) mean(FMRI{i_region,2}) mean(FMRI{i_region,3})],...
        [std(FMRI{i_region,1}) std(FMRI{i_region,2}) std(FMRI{i_region,3})]/2)
    errorbar([0.1 .5 1], [mean(FUSI{i_region,1}) mean(FUSI{i_region,2}) mean(FUSI{i_region,3})],...
        [std(FUSI{i_region,1}) std(FUSI{i_region,2}) std(FUSI{i_region,3})]/2,'r')
    title(T{i_region})
    set(gca,'xlim',[0 1.1]); %grid minor
    
    figure(26); subplot(2,2,i_region); hold on
    errorbar([0.1 .5 1], [mean(FMRI_noise{i_region,1}) mean(FMRI_noise{i_region,2}) mean(FMRI_noise{i_region,3})],...
        [std(FMRI_noise{i_region,1}) std(FMRI_noise{i_region,2}) std(FMRI_noise{i_region,3})]/2)
    errorbar([0.1 .5 1], [mean(FUSI_noise{i_region,1}) mean(FUSI_noise{i_region,2}) mean(FUSI_noise{i_region,3})],...
        [std(FUSI_noise{i_region,1}) std(FUSI_noise{i_region,2}) std(FUSI_noise{i_region,3})]/2,'r')
    title(T{i_region})
    set(gca,'xlim',[0 1.1]); %grid minor
    
    Zscore(1,i_region + 5) = {T{i_region}};
    Zscore(2:13,i_region + 5) = num2cell(horzcat(FMRI{i_region,:})');
    Zscore(14:25,i_region + 5) = num2cell(horzcat(FUSI{i_region,:})');
    
end

%%
% Write stat file for R - Z score fmri/fusi
Tzscore = ['C:\Users\Brad\Dropbox\fUSI_fMRI_Zscore.txt'];
fid = fopen(Tzscore,'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s\n', Zscore{1,:});
for K = 2:size(Zscore,1)
    fprintf(fid, '%.0f %s %.2f %.0f %.0f %.2f %.2f %.2f %.2f\n', Zscore{K,:});
end
fclose(fid)






%%
lfp_amp_file = [Dir 'Data_Processed' slash 'ephys' slash 'Thy1_LFP_AMP_values_GROUP.mat'];
load(lfp_amp_file); LFP = LFP_region_GRP;

Data = cell(13,6,3);
Data(1,:,1) = {'Animal ID','Intensity','R1','R2','R3','R4'};
Data(1,:,2) = {'Animal ID','Intensity','R1','R2','R3','R4'};
Data(1,:,3) = {'Animal ID','Intensity','R1','R2','R3','R4'};
% Data(1,:,4) = {'Animal ID','Intensity','R1','R2','R3','R4'};
% Data(1,:,5) = {'Animal ID','Intensity','R1','R2','R3','R4'};
Data(2:13,1,1) = {'1';'1';'1';'2';'2';'2';'3';'3';'3';'4';'4';'4'};
Data(2:13,1,2) = {'1';'1';'1';'2';'2';'2';'3';'3';'3';'4';'4';'4'};
Data(2:13,1,3) = {'1';'1';'1';'2';'2';'2';'3';'3';'3';'4';'4';'4'};
% Data(2:13,1,4) = {'1';'1';'1';'2';'2';'2';'3';'3';'3';'4';'4';'4'};
% Data(2:13,1,5) = {'1';'1';'1';'2';'2';'2';'3';'3';'3';'4';'4';'4'};
Data(2:13,2,1) = {'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';};
Data(2:13,2,2) = {'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';};
Data(2:13,2,3) = {'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';};
% Data(2:13,2,4) = {'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';};
% Data(2:13,2,5) = {'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';'0.1';'0.5';'1.0';};

% LFP Data Organization
LFP = reshape(cell2mat(LFP),[4 4 3]);
for i_region = 1:4
    for i_mouse = 1:4
        region_data_tot(:,i_mouse,i_region) = squeeze(LFP(i_mouse,i_region,:));
    end
    region_data = region_data_tot(:,:,i_region);
    Data(2:13,2 + i_region,1) = num2cell(region_data(:));
end

% FUSI Data Organization
ridx = [1 5 15 14];
for i_region = 1:4
    region_data = squeeze(fusi(ridx(i_region),:,:))';
    Data(2:13,2 + i_region,2) = num2cell(region_data(:));
end

% FMRI Data Organization
ridx = [1 5 15 14];
for i_region = 1:4
    region_data = squeeze(fmri(ridx(i_region),:,:))';
    Data(2:13,2 + i_region,3) = num2cell(region_data(:));
end


% FUSI Vascular Data Organization
% ridx = [1 5 15 14];
% for i_region = 1:4
%     region_data = squeeze(FUSIv(ridx(i_region),:,:))';
%     Data(2:13,2 + i_region,4) = num2cell(region_data(:));
%     region_data = squeeze(FUSIa(ridx(i_region),:,:))';
%     Data(2:13,2 + i_region,5) = num2cell(region_data(:));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT
% Ephys alone
figure; subplot(1,2,1); hold on
errorbar([0.1 0.5 1.0], mean(squeeze(LFP(2,:,:)),1), std(squeeze(LFP(2,:,:)),1)/2)
errorbar([0.1 0.5 1.0], mean(squeeze(LFP(4,:,:)),1), std(squeeze(LFP(4,:,:)),1)/2)
set(gca,'xlim',[0 1.1]); legend('L cortex','R cortex')
subplot(1,2,2); hold on
errorbar([0.1 0.5 1.0], mean(squeeze(LFP(1,:,:)),1), std(squeeze(LFP(1,:,:)),1)/2)
errorbar([0.1 0.5 1.0], mean(squeeze(LFP(3,:,:)),1), std(squeeze(LFP(3,:,:)),1)/2)
set(gca,'xlim',[0 1.1]); legend('L CPU','R CPU')

% FUSI
figure;
t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
for i_region = 1:4
    subplot(2,2,i_region); hold on
    Data1(:,i_region) = vertcat(Data{2:end,2+i_region,1});
    Data2(:,i_region) = vertcat(Data{2:end,2+i_region,2});
    scatter(Data1(1:3,i_region),Data2(1:3,i_region),'b','filled')
    scatter(Data1(4:6,i_region),Data2(4:6,i_region),'r','filled')
    scatter(Data1(7:9,i_region),Data2(7:9,i_region),'g','filled')
    scatter(Data1(10:12,i_region),Data2(10:12,i_region),'k','filled')

    mdl = fitlm(Data1(:,i_region),Data2(:,i_region));
    B = mdl.Coefficients.Estimate;                      % Coefficients
    CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
    [Ypred,YCI] = predict(mdl, Data1(:,i_region));      % Fitted Regression Line & Confidence Intervals

    plot(Data1(:,i_region), Ypred,'-r',Data1(:,i_region), YCI, '--r')
    title([t_bank{i_region} ': Rsq = ' num2str(mdl.Rsquared.Ordinary,3) ', p = ' num2str(mdl.Coefficients.pValue(2),3)]);
    if ismember(i_region,[1 3])
        ylim([-5 20]); xlim([0 125]);
    elseif ismember(i_region,[2 4])
        ylim([-20 150]); xlim([0 125]);
    end
end
supertitle('fUSI')

% FMRI
figure;
t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
for i_region = 1:4
    subplot(2,2,i_region); hold on
    Data1(:,i_region) = vertcat(Data{2:end,2+i_region,1});
    Data2(:,i_region) = vertcat(Data{2:end,2+i_region,3});
    scatter(Data1(1:3,i_region),Data2(1:3,i_region),'b','filled')
    scatter(Data1(4:6,i_region),Data2(4:6,i_region),'r','filled')
    scatter(Data1(7:9,i_region),Data2(7:9,i_region),'g','filled')
    scatter(Data1(10:12,i_region),Data2(10:12,i_region),'k','filled')

    mdl = fitlm(Data1(:,i_region),Data2(:,i_region));
    B = mdl.Coefficients.Estimate;                      % Coefficients
    CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
    [Ypred,YCI] = predict(mdl, Data1(:,i_region));      % Fitted Regression Line & Confidence Intervals

    plot(Data1(:,i_region), Ypred,'-r',Data1(:,i_region), YCI, '--r')
    title([t_bank{i_region} ': Rsq = ' num2str(mdl.Rsquared.Ordinary,3) ', p = ' num2str(mdl.Coefficients.pValue(2),3)]);
end
supertitle('fMRI')

% COMBINED PLOT
t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
for i_region = 1:4
    figure; subplot(2,2,1); hold on
    Data1(:,i_region) = vertcat(Data{2:end,2+i_region,1}); % LFP
    Data2(:,i_region) = vertcat(Data{2:end,2+i_region,2}); % fUSI
    Data3(:,i_region) = vertcat(Data{2:end,2+i_region,3}); % fMRI

    plot3(Data1(1:3,i_region),Data2(1:3,i_region),Data3(1:3,i_region),'linewidth',2,'color',[0 153 153]/255)
    plot3(Data1(4:6,i_region),Data2(4:6,i_region),Data3(4:6,i_region),'linewidth',2,'color',[150 0 0]/255)
    plot3(Data1(7:9,i_region),Data2(7:9,i_region),Data3(7:9,i_region),'linewidth',2,'color',[204 102 0]/255)
    plot3(Data1(10:12,i_region),Data2(10:12,i_region),Data3(10:12,i_region),'linewidth',2,'color',[32 32 32]/255)

    scatter3(Data1(1:3:end,i_region),Data2(1:3:end,i_region),Data3(1:3:end,i_region),150,...
        [[0 153 153]/255;[150 0 0]/255;[204 102 0]/255;[32 32 32]/255],'filled','marker','o')
    scatter3(Data1(2:3:end,i_region),Data2(2:3:end,i_region),Data3(2:3:end,i_region),150,...
        [[0 153 153]/255;[150 0 0]/255;[204 102 0]/255;[32 32 32]/255],'filled','marker','s')
    scatter3(Data1(3:3:end,i_region),Data2(3:3:end,i_region),Data3(3:3:end,i_region),150,...
        [[0 153 153]/255;[150 0 0]/255;[204 102 0]/255;[32 32 32]/255],'filled','marker','d')


    xlabel('LFP'); ylabel('fUSI'); zlabel('fMRI')

    if ismember(i_region,[1 3])
        xmin = -25; xmax = 15000;
        ymin = -25; ymax = 35; ymin2 = 0; ymax2 = 20;
        zmin = -.5; zmax = 2; zmin2 = 0; zmax2 = 1;
    elseif ismember(i_region,[2 4])
        xmin = -25; xmax = 15000;
        ymin = -50; ymax = 235; ymin2 = 0; ymax2 = 210;
        zmin = -1; zmax = 2.5; zmin2 = 0; zmax2 = 2;
    end
    xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]);

    % LFP vs fUSI
    S = 75; L = 1.25;
    plot3(Data1(1:3,i_region),Data2(1:3,i_region),zmin*ones(3,1),'linewidth',L,'color',[102 255 255]/255)
    plot3(Data1(4:6,i_region),Data2(4:6,i_region),zmin*ones(3,1),'linewidth',L,'color',[255 102 102]/255)
    plot3(Data1(7:9,i_region),Data2(7:9,i_region),zmin*ones(3,1),'linewidth',L,'color',[255 178 102]/255)
    plot3(Data1(10:12,i_region),Data2(10:12,i_region),zmin*ones(3,1),'linewidth',L,'color',[128 128 128]/255)

    scatter3(Data1(1:3:end,i_region),Data2(1:3:end,i_region),zmin*ones(4,1),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter3(Data1(2:3:end,i_region),Data2(2:3:end,i_region),zmin*ones(4,1),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter3(Data1(3:3:end,i_region),Data2(3:3:end,i_region),zmin*ones(4,1),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')

	% LFP vs fMRI
    plot3(Data1(1:3,i_region),ymax*ones(3,1),Data3(1:3,i_region),'linewidth',L,'color',[102 255 255]/255)
    plot3(Data1(4:6,i_region),ymax*ones(3,1),Data3(4:6,i_region),'linewidth',L,'color',[255 102 102]/255)
    plot3(Data1(7:9,i_region),ymax*ones(3,1),Data3(7:9,i_region),'linewidth',L,'color',[255 178 102]/255)
    plot3(Data1(10:12,i_region),ymax*ones(3,1),Data3(10:12,i_region),'linewidth',L,'color',[128 128 128]/255)

    scatter3(Data1(1:3:end,i_region),ymax*ones(4,1),Data3(1:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter3(Data1(2:3:end,i_region),ymax*ones(4,1),Data3(2:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter3(Data1(3:3:end,i_region),ymax*ones(4,1),Data3(3:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')

	% fUSI vs fMRI
    plot3(xmin*ones(3,1),Data2(1:3,i_region),Data3(1:3,i_region),'linewidth',L,'color',[102 255 255]/255)
    plot3(xmin*ones(3,1),Data2(4:6,i_region),Data3(4:6,i_region),'linewidth',L,'color',[255 102 102]/255)
    plot3(xmin*ones(3,1),Data2(7:9,i_region),Data3(7:9,i_region),'linewidth',L,'color',[255 178 102]/255)
    plot3(xmin*ones(3,1),Data2(10:12,i_region),Data3(10:12,i_region),'linewidth',L,'color',[128 128 128]/255)

    scatter3(xmin*ones(4,1),Data2(1:3:end,i_region),Data3(1:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter3(xmin*ones(4,1),Data2(2:3:end,i_region),Data3(2:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter3(xmin*ones(4,1),Data2(3:3:end,i_region),Data3(3:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')

    grid on
    ax = gca;
    set(ax,'GridColor',[0 0 0],'linewidth',1,'color',[1 1 1])
%     set(gcf,'position',[2175 41 790 740]);
    view(38,35);
    
    % Individual plots
    m11 = mean(Data1(1:3:end,i_region)); m12 = mean(Data1(2:3:end,i_region)); m13 = mean(Data1(3:3:end,i_region)); % ephys
    m21 = mean(Data2(1:3:end,i_region)); m22 = mean(Data2(2:3:end,i_region)); m23 = mean(Data2(3:3:end,i_region)); % fusi
    m31 = mean(Data3(1:3:end,i_region)); m32 = mean(Data3(2:3:end,i_region)); m33 = mean(Data3(3:3:end,i_region)); % fmri

    s11 = std(Data1(1:3:end,i_region))/2; s12 = std(Data1(2:3:end,i_region))/2; s13 = std(Data1(3:3:end,i_region))/2;
    s21 = std(Data2(1:3:end,i_region))/2; s22 = std(Data2(2:3:end,i_region))/2; s23 = std(Data2(3:3:end,i_region))/2;
    s31 = std(Data3(1:3:end,i_region))/2; s32 = std(Data3(2:3:end,i_region))/2; s33 = std(Data3(3:3:end,i_region))/2;
    
    
    subplot(2,2,2); hold on
    % LFP vs fUSI
    S = 75; L = 1.25;
    plot(Data1(1:3,i_region),Data2(1:3,i_region),'linewidth',L,'color',[102 255 255]/255)
    plot(Data1(4:6,i_region),Data2(4:6,i_region),'linewidth',L,'color',[255 102 102]/255)
    plot(Data1(7:9,i_region),Data2(7:9,i_region),'linewidth',L,'color',[255 178 102]/255)
    plot(Data1(10:12,i_region),Data2(10:12,i_region),'linewidth',L,'color',[128 128 128]/255)
    % average
    plot([m11 m12 m13],[m21 m22 m23],'linewidth',L,'color','k')
    
    scatter(Data1(1:3:end,i_region),Data2(1:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter(Data1(2:3:end,i_region),Data2(2:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter(Data1(3:3:end,i_region),Data2(3:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')
    % average
    scatter([m11 m12 m13],[m21 m22 m23],'k','filled','marker','s')
    
    xlabel('LFP'); ylabel('fUSI'); grid on; xlim([xmin xmax]); ylim([ymin2 ymax2])
    
    % LFP vs fMRI
    subplot(2,2,3); hold on
    plot(Data1(1:3,i_region),Data3(1:3,i_region),'linewidth',L,'color',[102 255 255]/255)
    plot(Data1(4:6,i_region),Data3(4:6,i_region),'linewidth',L,'color',[255 102 102]/255)
    plot(Data1(7:9,i_region),Data3(7:9,i_region),'linewidth',L,'color',[255 178 102]/255)
    plot(Data1(10:12,i_region),Data3(10:12,i_region),'linewidth',L,'color',[128 128 128]/255)
    % average
    plot([m11 m12 m13],[m31 m32 m33],'linewidth',L,'color','k')
    
    scatter(Data1(1:3:end,i_region),Data3(1:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter(Data1(2:3:end,i_region),Data3(2:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter(Data1(3:3:end,i_region),Data3(3:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')
    % average
    scatter([m11 m12 m13],[m31 m32 m33],'k','filled','marker','s')
    
    xlabel('LFP'); ylabel('fMRI'); grid on; xlim([xmin xmax]); ylim([zmin2 zmax2])
    
    % fUSI vs fMRI
    subplot(2,2,4); hold on
    plot(Data2(1:3,i_region),Data3(1:3,i_region),'linewidth',L,'color',[102 255 255]/255)
    plot(Data2(4:6,i_region),Data3(4:6,i_region),'linewidth',L,'color',[255 102 102]/255)
    plot(Data2(7:9,i_region),Data3(7:9,i_region),'linewidth',L,'color',[255 178 102]/255)
    plot(Data2(10:12,i_region),Data3(10:12,i_region),'linewidth',L,'color',[128 128 128]/255)
    % average
    plot([m21 m22 m23],[m31 m32 m33],'linewidth',L,'color','k')
    
    scatter(Data2(1:3:end,i_region),Data3(1:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter(Data2(2:3:end,i_region),Data3(2:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter(Data2(3:3:end,i_region),Data3(3:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')
    % average
    scatter([m21 m22 m23],[m31 m32 m33],'k','filled','marker','s')
    xlabel('fUSI'); ylabel('fMRI'); grid on
    supertitle(t_bank{i_region})
    
    % Plotting only average for LFP vs fUSI/fMRI
    figure; hold on
    xlabel('LFP'); grid on
    % fMRI
    yyaxis left; ylabel('fMRI'); 
    fmriC = [179 205 227; 140 150 198; 136 65 157]/255;
%     plot([m11 m12 m13],[m31 m32 m33],'linewidth',L,'color','k')
    errorbar([m11 m12 m13],[m31 m32 m33],[s31 s32 s33],[s31 s32 s33],[s11 s12 s13],[s11 s12 s13],'k')
    scatter([m11 m12 m13],[m31 m32 m33],150,fmriC,'filled','marker','s')
    if ismember(i_region,[1 3])
        ylim([0 .5])
    elseif ismember(i_region,[2,4])
        ylim([0 2])
    end
    % fUSI
    yyaxis right; ylabel('fUSI'); 
    fusiC = [175 216 166; 90 184 70; 32 104 52]/255;
%     plot([m11 m12 m13],[m21 m22 m23],'linewidth',L,'color','b')
    errorbar([m11 m12 m13],[m21 m22 m23],[s21 s22 s23],[s21 s22 s23],[s11 s12 s13],[s11 s12 s13],'b')
    scatter([m11 m12 m13],[m21 m22 m23],150,fusiC,'filled','marker','s')
    if ismember(i_region,[1 3])
        ylim([0 20])
    elseif ismember(i_region,[2,4])
        ylim([0 200])
    end
    title(t_bank{i_region})
    
    if ismember(i_region,[1:4])
        fprintf('\nRegion: %s\n,',t_bank{i_region})
        [m21 m22 m23]./[m31 m32 m33]
        [m21 m22 m23]
        [m31 m32 m33]
        
    end
end

%%
% FUSI Vascular Components
% COMBINED PLOT
t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
for i_region = 1:4
    figure; hold on
    Data1(:,i_region) = vertcat(Data{2:end,2+i_region,4}); % fUSI Ven
    Data2(:,i_region) = vertcat(Data{2:end,2+i_region,5}); % fUSI Art
    Data3(:,i_region) = vertcat(Data{2:end,2+i_region,1}); % LFP 

    plot3(Data1(1:3,i_region),Data2(1:3,i_region),Data3(1:3,i_region),'linewidth',2,'color',[0 153 153]/255)
    plot3(Data1(4:6,i_region),Data2(4:6,i_region),Data3(4:6,i_region),'linewidth',2,'color',[150 0 0]/255)
    plot3(Data1(7:9,i_region),Data2(7:9,i_region),Data3(7:9,i_region),'linewidth',2,'color',[204 102 0]/255)
    plot3(Data1(10:12,i_region),Data2(10:12,i_region),Data3(10:12,i_region),'linewidth',2,'color',[32 32 32]/255)

    scatter3(Data1(1:3:end,i_region),Data2(1:3:end,i_region),Data3(1:3:end,i_region),150,...
        [[0 153 153]/255;[150 0 0]/255;[204 102 0]/255;[32 32 32]/255],'filled','marker','o')
    scatter3(Data1(2:3:end,i_region),Data2(2:3:end,i_region),Data3(2:3:end,i_region),150,...
        [[0 153 153]/255;[150 0 0]/255;[204 102 0]/255;[32 32 32]/255],'filled','marker','s')
    scatter3(Data1(3:3:end,i_region),Data2(3:3:end,i_region),Data3(3:3:end,i_region),150,...
        [[0 153 153]/255;[150 0 0]/255;[204 102 0]/255;[32 32 32]/255],'filled','marker','d')


    xlabel('fUSI Vein'); ylabel('fUSI Artery'); zlabel('LFP')

    if ismember(i_region,[1 3])
        xmin = -25; xmax = 35;
        ymin = -25; ymax = 35;
        zmin = -25; zmax = 150;
    elseif ismember(i_region,[2 4])
        xmin = -50; xmax = 175;
        ymin = -50; ymax = 175;
        zmin = -25; zmax = 150;
    end
    xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]);

    % LFP vs fUSI
    S = 75; L = 1.25;
    plot3(Data1(1:3,i_region),Data2(1:3,i_region),zmin*ones(3,1),'linewidth',L,'color',[102 255 255]/255)
    plot3(Data1(4:6,i_region),Data2(4:6,i_region),zmin*ones(3,1),'linewidth',L,'color',[255 102 102]/255)
    plot3(Data1(7:9,i_region),Data2(7:9,i_region),zmin*ones(3,1),'linewidth',L,'color',[255 178 102]/255)
    plot3(Data1(10:12,i_region),Data2(10:12,i_region),zmin*ones(3,1),'linewidth',L,'color',[128 128 128]/255)

    scatter3(Data1(1:3:end,i_region),Data2(1:3:end,i_region),zmin*ones(4,1),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter3(Data1(2:3:end,i_region),Data2(2:3:end,i_region),zmin*ones(4,1),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter3(Data1(3:3:end,i_region),Data2(3:3:end,i_region),zmin*ones(4,1),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')

    % LFP vs fMRI
    plot3(Data1(1:3,i_region),ymax*ones(3,1),Data3(1:3,i_region),'linewidth',L,'color',[102 255 255]/255)
    plot3(Data1(4:6,i_region),ymax*ones(3,1),Data3(4:6,i_region),'linewidth',L,'color',[255 102 102]/255)
    plot3(Data1(7:9,i_region),ymax*ones(3,1),Data3(7:9,i_region),'linewidth',L,'color',[255 178 102]/255)
    plot3(Data1(10:12,i_region),ymax*ones(3,1),Data3(10:12,i_region),'linewidth',L,'color',[128 128 128]/255)

    scatter3(Data1(1:3:end,i_region),ymax*ones(4,1),Data3(1:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter3(Data1(2:3:end,i_region),ymax*ones(4,1),Data3(2:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter3(Data1(3:3:end,i_region),ymax*ones(4,1),Data3(3:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')

    % fUSI vs fMRI
    plot3(xmin*ones(3,1),Data2(1:3,i_region),Data3(1:3,i_region),'linewidth',L,'color',[102 255 255]/255)
    plot3(xmin*ones(3,1),Data2(4:6,i_region),Data3(4:6,i_region),'linewidth',L,'color',[255 102 102]/255)
    plot3(xmin*ones(3,1),Data2(7:9,i_region),Data3(7:9,i_region),'linewidth',L,'color',[255 178 102]/255)
    plot3(xmin*ones(3,1),Data2(10:12,i_region),Data3(10:12,i_region),'linewidth',L,'color',[128 128 128]/255)

    scatter3(xmin*ones(4,1),Data2(1:3:end,i_region),Data3(1:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','o')
    scatter3(xmin*ones(4,1),Data2(2:3:end,i_region),Data3(2:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','s')
    scatter3(xmin*ones(4,1),Data2(3:3:end,i_region),Data3(3:3:end,i_region),S,...
        [[102 255 255]/255;[255 102 102]/255;[255 178 102]/255;[128 128 128]/255],'filled','marker','d')

    grid on
    ax = gca;
    set(ax,'GridColor',[0 0 0],'linewidth',1,'color',[1 1 1])
    set(gcf,'position',[2175 41 790 740]);

    view(38,35);
    title(t_bank{i_region})

end





