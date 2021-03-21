function fmri_plot_timeseries_analysis(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
stimC = [179 205 227; 140 150 198; 136 65 157]/255;

fmriroi = param.fmriroi;
descr = param.descr;

h1 = figure(1); clf
h2 = figure(2); clf

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
            if strcmp(fmriroi{i_roi},'RM1')
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
                tmp = tmp - repmat(mean(ts,2),[size(tmp,1) size(tmp,2)]);
                tmp = tmp - repmat(mean(ts1(1:10)),[size(tmp,1) size(tmp,2)]); 
                
                stdshade(tmp,0.25,stimC(i_stim,:))
%                 plot(ts,'color',stimC(i_stim,:),'linewidth',2)
            end
            
             % Plot peak time series voxel-wise for ROIs
            if contains(['RM1' 'LM1'] ,fmriroi{i_roi})
                tmp = nanmean(cell2mat(tsnormaveZ),2);
                tmp = reshape(tmp,[size(grp_data,1) size(grp_data,2)]);
                tmp2 = zeros(size(grp_data,1),size(grp_data,2));
                tmp2(pts{i_roi,1}) = tmp(pts{i_roi,1});
                tmp2 = imresize(tmp2,size(anat));
                anat(tmp2 > .05) = tmp2 (tmp2 > .05);
                figure(i_stim+250); imagesc(anat); colormap(Cbar); caxis(CAX)
            end
            
            
        end
        
    end
end
figure(h1); supertitle('ROI time series: Intensity')
saveas(h1,[storage descr '_ROI_TS_fMRI.svg']);
saveas(h1,[storage descr '_ROI_TS_fMRI.fig']);

figure(h2); title('RM1 Full Timeseries')
saveas(h2,[storage descr '_RM1_Full_Timeseries_fMRI.svg']);
saveas(h2,[storage descr '_RM1_Full_Timeseries_fMRI.fig']);
