function fmri_plot_timeseries_analysis(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
stimC = [179 205 227; 140 150 198; 136 65 157]/255;

fmriroi = param.fmriroi;
descr = param.descr;

h1 = figure(1); clf
h2 = figure(2); clf

SingleROI = 'RCPu';

for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_time_series_info.mat'];
    grp_file = [stim_storage stim{i_stim} '_grp_data.mat'];
    count_file = [stim_storage stim{i_stim} '_active_voxel_count.mat'];
    
    [anat,Cbar,CAX] = fus_prep_anat(storage,slash,param,-100,0,20);
    if exist(stim_file,'file') && exist(count_file,'file') && exist(grp_file,'file')
        
        load(stim_file)
        load(count_file)
        load(grp_file)
        for i_roi = 1:20
            
            % Anatomical ROI timeseries (across intensities)
            figure(h1); subplot(5,4,i_roi); hold on
            stdshade(vertcat(tsnormave_roi{i_roi,:}),0.25,stimC(i_stim,:))
            y1=get(gca,'ylim');
            plot([11 11],y1,'k'); plot([19 19],y1,'k')
            title(fmriroi{i_roi})
            
            % Full timeseries for local stimulation site
            if strcmp(fmriroi{i_roi},SingleROI)
                figure(h2); hold on
                if isequal(i_stim,1)
                    start = [30 70 110 150 190];
                    for i = 1:size(start,2)
                        x = [start(i) start(i)+8 start(i)+8 start(i)];
                        y = [-1 -1 2 2];
                        patch(x,y,'k','facealpha',.25)
                    end
                end
                ts = mean(vertcat(tsnorm_roi{i_roi,:}),1);
                ts1 = ts - mean(ts,2);
                ts2 = ts1 - mean(ts1(1:10));
                
                tmp = vertcat(tsnorm_roi{i_roi,:});
                for i= 1:size(tmp,1)
                    tmp(1,:) = detrend(tmp(1,:));
                end
                tmp = tmp - repmat(mean(mean(tmp,2),1),[size(tmp,1) size(tmp,2)]);
                tmp = tmp - repmat(mean(mean(tmp(1:10)),1),[size(tmp,1) size(tmp,2)]); 
                
                stdshade(tmp,0.25,stimC(i_stim,:))
                set(gca,'ylim',[-2 2])
%                 plot(ts,'color',stimC(i_stim,:),'linewidth',2)
            end
            
             % Plot peak time series voxel-wise for ROIs
            if contains(['RM1' 'LM1'] ,fmriroi{i_roi})
                tmp = nanmean(cell2mat(tsnormavepeak),2);
                tmp = reshape(tmp,[size(grp_data,1) size(grp_data,2)]);
                tmp2 = zeros(size(grp_data,1),size(grp_data,2));
                tmp2(pts{i_roi,1}) = tmp(pts{i_roi,1});
                tmp2 = imresize(tmp2,size(anat));
                anat(tmp2 > .05) = tmp2 (tmp2 > .05);
                figure(i_stim+250); imagesc(anat); colormap(Cbar); caxis(CAX)
            end
            
            
        end
        
    end
    
    % Plot noise z-score
    NOISE_Z{i_stim} = cell2mat(tsnormaveZ_noise(end,:));
    ROI_Z{i_stim} = tsnormaveZ_roi;
end
figure(h1); supertitle('ROI time series: Intensity')
saveas(h1,[storage descr '_ROI_TS_fMRI.svg']);
saveas(h1,[storage descr '_ROI_TS_fMRI.fig']);

figure(h2); title('RM1 Full Timeseries')
saveas(h2,[storage descr '_' SingleROI '_Full_Timeseries_fMRI.svg']);
saveas(h2,[storage descr '_' SingleROI '_Full_Timeseries_fMRI.fig']);


%% ROI CIRCUIT Z SCORES

roi_int = [];
roi_int = [roi_int find(contains(fmriroi, 'RM1'))];
roi_int = [roi_int find(contains(fmriroi, 'LM1'))];
roi_int = [roi_int find(contains(fmriroi, 'RCPu'))];
roi_int = [roi_int find(contains(fmriroi, 'LCPu'))];

for i = 1:3
    ROIZ(:,:,i) = cell2mat(ROI_Z{i});
end


f = figure(62);
num_an = size(tsnorm_roi,2);
for i_roi = 1:size(roi_int,2)
    subplot(2,2,i_roi)
    % Plot Z score for ROI
    D = ROIZ(roi_int(i_roi),:,:);
    D = reshape(D,[],1);
    Dl = [repmat('1',num_an,1);...
        repmat('2',num_an,1);...
        repmat('3',num_an,1)];
    boxplot(D,Dl);
    set(gca,'xticklabel',{'0.1' '0.5' '1.0'},'ylim',[-10 35])
    title(fmriroi{roi_int(i_roi)})
end
saveas(f,[storage descr '_ROI_CIRCUIT_Z.fig']);
saveas(f,[storage descr '_ROI_CIRCUIT_Z.svg']);

%% NOISE Z SCORES
% PLOT NOISE ROI Z-SCORE INFO
h27 = figure;
subplot(2,1,1);
D = [NOISE_Z{1}' NOISE_Z{2}' NOISE_Z{3}'];
Dl = [repmat('1',size(D,1),1);...
    repmat('2',size(D,1),1);...
    repmat('3',size(D,1),1)];
boxplot(D,Dl);
set(gca,'xticklabel',{'0.1' '0.5' '1.0'})
hold on
Xr = ones(size(D)).*(1+(rand(size(D))-0.5)/10);
Xr(:,2) = Xr(:,2) + 1; Xr(:,3) = Xr(:,3) + 2;
scatter(Xr(:),D(:),'r','filled')
title('raw noise Z')

subplot(2,1,2);
boxplot(abs(D),Dl);
hold on
Xr = ones(size(D)).*(1+(rand(size(D))-0.5)/10);
Xr(:,2) = Xr(:,2) + 1; Xr(:,3) = Xr(:,3) + 2;
scatter(Xr(:),abs(D(:)),'r','filled')
set(gca,'xticklabel',{'0.1' '0.5' '1.0'},'ylim',[-1 10])
title('abs noise Z')
figure(h27); supertitle('Noise Z score')
saveas(h27,[storage descr '_NOISE_Z.fig']);
saveas(h27,[storage descr '_NOISE_Z.svg']);
