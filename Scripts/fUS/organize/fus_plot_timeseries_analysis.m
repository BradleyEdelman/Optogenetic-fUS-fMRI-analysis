function fus_plot_timeseries_analysis(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
stimC = [175 216 166; 90 184 70; 32 104 52]/255;
redC = [252 174 145; 222 45 38; 165 15 21]/255;
blueC = [189 215 231; 49 130 189; 8 81 156]/255;
vascC = [1 0 0; 0 0 1];

fusroi = param.fusroi;
descr = param.descr;

h1 = figure(1); clf
h2 = figure(2); clf
h3 = figure(3); clf
h4 = figure(4); clf
h5 = figure(5); clf
h6 = figure(6); clf
h7 = figure(7); clf
h8 = figure(8); clf

StatAUC(1,1:5) = {'subjid','vessel','intensity','pairInt','pairVes'};
id = repmat(1:size(base_fold,1),3,1); id = repmat(id(:),2,1);
StatAUC(2:1+size(base_fold,1)*6,1) = num2cell(id);
StatAUC(2:1+size(base_fold,1)*3,2) = {'Art'};
StatAUC(2+size(base_fold,1)*3:1+size(base_fold,1)*6,2) = {'Vein'};
StatAUC(2:1+size(base_fold,1)*6,3) = num2cell(repmat([1;2;3],2*size(base_fold,1),1));
StatAUC(2:1+size(base_fold,1)*6,4) = num2cell(repmat([1;2;3],2*size(base_fold,1),1));
StatAUC(2:1+size(base_fold,1)*3,5) = num2cell(repmat([1;2;3],size(base_fold,1),1));
StatAUC(2+size(base_fold,1)*3:1+size(base_fold,1)*6,5) = num2cell(repmat([4;5;6],size(base_fold,1),1));
StatPk = StatAUC; StatCt = StatAUC;

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
            title(fusroi{i_roi})
            
            % Full timeseries for local stimulation site
            if strcmp(fusroi{i_roi},'RM1')
                figure(h8); hold on
                if isequal(i_stim,3)
                    start = [30 69 108 147 186];
                    for i = 1:size(start,2)
                        x = [start(i) start(i)+8 start(i)+8 start(i)];
                        y = [-25 -25 100 100];
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
%                 plot(ts2,'color',stimC(i_stim,:),'linewidth',2)
            end
            
            % Plot peak time series voxel-wise for ROIs
            if contains(['RM1' 'LM1' 'RCPu' 'LCPu'] ,fusroi{i_roi})
                tmp = nanmean(cell2mat(tsnormaveZ),2);
                tmp = reshape(tmp,[size(grp_data,1) size(grp_data,2)]);
                tmp2 = zeros(size(grp_data,1),size(grp_data,2));
                tmp2(pts{i_roi,1}) = tmp(pts{i_roi,1});
                tmp2 = imresize(tmp2,size(anat));
                anat(tmp2 > 0.05) = tmp2 (tmp2 > 0.05);
                figure(i_stim+250); imagesc(anat); colormap(Cbar); caxis(CAX)
            end
            
            % Sig active arterial ROI timeseries (across intensities)
            figure(h2); subplot(5,4,i_roi); hold on
            stdshade(vertcat(tsnormave_roi_sigpos{i_roi,:}),0.25,redC(i_stim,:))
            y1=get(gca,'ylim');
            plot([11 11],y1,'k'); plot([19 19],y1,'k')
            title(fusroi{i_roi})

            % Sig active venous ROI timeseries (across intensities)
            figure(h3); subplot(5,4,i_roi); hold on
            stdshade(vertcat(tsnormave_roi_signeg{i_roi,:}),0.25,blueC(i_stim,:))
            y1=get(gca,'ylim');
            plot([11 11],y1,'k'); plot([19 19],y1,'k')
            title(fusroi{i_roi})
            
            % Sig active arterial vs venous timeseries (within intensity)
            figure(h4); subplot(5,4,i_roi); cla; hold on
            stdshade(vertcat(tsnormave_roi_sigpos{i_roi,:}),0.25,vascC(1,:))
            stdshade(vertcat(tsnormave_roi_signeg{i_roi,:}),0.25,vascC(2,:))
            y1=get(gca,'ylim');
            plot([11 11],y1,'k'); plot([19 19],y1,'k')
            title(fusroi{i_roi})
            
            % AUC for arterial and venous components
            AUC_pos_mean{i_roi,i_stim} = nanmean(vertcat(tsnormaveAUC_roi_sigpos{i_roi,:}));
            AUC_pos_sem{i_roi,i_stim} = nanstd(vertcat(tsnormaveAUC_roi_sigpos{i_roi,:}))/...
                sqrt(sum(~isnan(vertcat(tsnormaveAUC_roi_sigpos{i_roi,:}))));
            AUC_neg_mean{i_roi,i_stim} = nanmean(vertcat(tsnormaveAUC_roi_signeg{i_roi,:}));
            AUC_neg_sem{i_roi,i_stim} = nanstd(vertcat(tsnormaveAUC_roi_signeg{i_roi,:}))/...
                sqrt(sum(~isnan(vertcat(tsnormaveAUC_roi_signeg{i_roi,:}))));
            
            % Format AUC values for stat testing (R .txt file)
            StatAUC{1,5+i_roi} = fusroi{i_roi};
            StatAUC(i_stim+1:3:1+size(base_fold,1)*3,5+i_roi) = num2cell(vertcat(tsnormaveAUC_roi_sigpos{i_roi,:}));
            StatAUC(i_stim+1+size(base_fold,1)*3:3:1+size(base_fold,1)*6,5+i_roi) = num2cell(vertcat(tsnormaveAUC_roi_signeg{i_roi,:}));

            % Peak activation for arterial and venous components
            peak_pos_mean{i_roi,i_stim} = nanmean(vertcat(tsnormavepeak_roi_sigpos{i_roi,:}));
            peak_pos_sem{i_roi,i_stim} = nanstd(vertcat(tsnormavepeak_roi_sigpos{i_roi,:}))/...
                sqrt(sum(~isnan(vertcat(tsnormavepeak_roi_sigpos{i_roi,:}))));
            peak_neg_mean{i_roi,i_stim} = nanmean(vertcat(tsnormavepeak_roi_signeg{i_roi,:}));
            peak_neg_sem{i_roi,i_stim} = nanstd(vertcat(tsnormavepeak_roi_signeg{i_roi,:}))/...
                sqrt(sum(~isnan(vertcat(tsnormavepeak_roi_signeg{i_roi,:}))));
            
            % Format peak act values for stat testing (R .txt file)
            StatPk{1,5+i_roi} = fusroi{i_roi};
            StatPk(i_stim+1:3:1+size(base_fold,1)*3,5+i_roi) = num2cell(vertcat(tsnormavepeak_roi_sigpos{i_roi,:}));
            StatPk(i_stim+1+size(base_fold,1)*3:3:1+size(base_fold,1)*6,5+i_roi) = num2cell(vertcat(tsnormavepeak_roi_signeg{i_roi,:}));
            
            % Pixel count for arterial and venous components
            clear sigpos_ct signeg_ct
            for i_mouse = 1:size(base_fold,1)
                sigpos_ct(i_mouse) = size(flow_roi_pos_sig_pts{i_roi,i_mouse},1);
                signeg_ct(i_mouse) = size(flow_roi_neg_sig_pts{i_roi,i_mouse},1);
            end
            ct_sigpos_mean{i_roi,i_stim} = nanmean(sigpos_ct);
            ct_sigpos_sem{i_roi,i_stim} = nanstd(sigpos_ct)/sum(size(~isnan(sigpos_ct),2));
            ct_signeg_mean{i_roi,i_stim} = nanmean(signeg_ct);
            ct_signeg_sem{i_roi,i_stim} = nanstd(signeg_ct)/sum(size(~isnan(signeg_ct),2));
            
            % Format pixel count values for stat testing (R .txt file)
            StatCt{1,5+i_roi} = fusroi{i_roi};
            StatCt(i_stim+1:3:1+size(base_fold,1)*3,5+i_roi) = num2cell(sigpos_ct');
            StatCt(i_stim+1+size(base_fold,1)*3:3:1+size(base_fold,1)*6,5+i_roi) = num2cell(signeg_ct');
            
        end
    end
    
    figure(h4); supertitle(['Arteriole vs Venule Time Series: ' stim{i_stim}]);
    saveas(h4,[stim_storage 'Art_vs_Ven_TS.svg']);
    saveas(h4,[stim_storage 'Art_vs_Ven_TS.fig']);
    
    
    TMP{i_stim,1} = tsnorm_peaks_pos;
    TMP{i_stim,2} = tsnorm_peaks_neg;
    
end
figure(h1); supertitle('ROI time series: Intensity')
saveas(h1,[storage descr '_ROI_TS.svg']);
saveas(h1,[storage descr '_ROI_TS.fig']);

figure(h2); supertitle('Arteriole time series: Instensity')
saveas(h2,[storage descr '_Artery_TS.svg']);
saveas(h2,[storage descr '_Artery_TS.fig']);

figure(h3); supertitle('Venule time series: Intensity')
saveas(h3,[storage descr '_Venule_TS.svg']);
saveas(h3,[storage descr '_Venule_TS.fig']);

if isunix
    Tvasauc = ['/home/bradley/Dropbox/fUSI_Vascular_AUC_' descr '.txt'];
    Tvaspk = ['/home/bradley/Dropbox/fUSI_Vascular_Pk_' descr '.txt'];
    Tvasct = ['/home/bradley/Dropbox/fUSI_Vascular_Ct_' descr '.txt'];
elseif ispc
    Tvasauc = ['C:\Users\Brad\Dropbox\fUSI_Vascular_AUC_' descr '.txt'];
    Tvaspk = ['C:\Users\Brad\Dropbox\fUSI_Vascular_Pk_' descr '.txt'];
    Tvasct = ['C:\Users\Brad\Dropbox\fUSI_Vascular_Ct_' descr '.txt'];
end

% Write stat file for R - AUC vascular components
StatAUC(:,[11,25]) = []; % Remove S2 regions!
fid = fopen(Tvasauc,'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatAUC{1,:});
for K = 2:size(StatAUC,1)
    fprintf(fid, '%.0f %s %.2f %.0f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatAUC{K,:});
end
fclose(fid)

% Write stat file for R - Peak Response vascular components
StatPk(:,[11,25]) = [];
fid = fopen(Tvaspk,'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatPk{1,:});
for K = 2:size(StatPk,1)
    fprintf(fid, '%.0f %s %.2f %.0f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatPk{K,:});
end
fclose(fid)

% Write stat file for R - Vascular active pixel count
StatCt(:,[11,25]) = [];
fid = fopen(Tvasct,'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatCt{1,:});
for K = 2:size(StatCt,1)
    fprintf(fid, '%.0f %s %.2f %.0f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatCt{K,:});
end
fclose(fid)

for i_roi = 1:20
    
    % Plot AUC for arterial and venous components
    figure(h5); subplot(5,4,i_roi); cla; hold on
    Ave = [vertcat(AUC_pos_mean{i_roi,:}) vertcat(AUC_neg_mean{i_roi,:})];
    Sem = [vertcat(AUC_pos_sem{i_roi,:}) vertcat(AUC_neg_sem{i_roi,:})];
    Sem(isnan(Sem)) = 0;
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
    title(fusroi{i_roi})
    
    % Plot peak activation for arterial and venous components
    figure(h6); subplot(5,4,i_roi); cla; hold on
    Ave = [vertcat(peak_pos_mean{i_roi,:}) vertcat(peak_neg_mean{i_roi,:})];
    Sem = [vertcat(peak_pos_sem{i_roi,:}) vertcat(peak_neg_sem{i_roi,:})];
    Sem(isnan(Sem)) = 0;
    ctrs = 1:2; hBar{1} = bar(ctrs, Ave'); drawnow; clear ctr ydt
    for k1 = 1:size(Ave',2)
        ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
        ydt(k1,:) = hBar{1}(k1).YData;
    end
    errorbar(ctr, ydt, zeros(2,3)', Sem, zeros(2,3)', zeros(2,3)', '.k',...
        'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
    hBar{2} = bar(ctrs, Ave'); delete(hBar{1})
    if max(Ave(:)) < 50; M = 100; else M = 200; end
    set(gca,'ylim',[0 M],'xtick',[1 2],'xticklabel',{'art' 'ven'});
    title(fusroi{i_roi})
    
    % Plot pixel count for arterial and venous components
    figure(h7); subplot(5,4,i_roi); cla; hold on
    Ave = [vertcat(ct_sigpos_mean{i_roi,:}) vertcat(ct_signeg_mean{i_roi,:})];
    Sem = [vertcat(ct_sigpos_sem{i_roi,:}) vertcat(ct_signeg_sem{i_roi,:})];
    Sem(isnan(Sem)) = 0;
    ctrs = 1:2; hBar{1} = bar(ctrs, Ave'); drawnow; clear ctr ydt
    for k1 = 1:size(Ave',2)
        ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
        ydt(k1,:) = hBar{1}(k1).YData;
    end
    errorbar(ctr, ydt, zeros(2,3)', Sem, zeros(2,3)', zeros(2,3)', '.k',...
        'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
    hBar{2} = bar(ctrs, Ave'); delete(hBar{1})
    if max(Ave(:)) < 50; M = 75; else M = 300; end
    set(gca,'ylim',[0 M],'xtick',[1 2],'xticklabel',{'art' 'ven'});
    title(fusroi{i_roi})
    
end
figure(h5); supertitle('Arteriole vs Venule AUC: Intensity')
saveas(h5,[storage descr '_Art_vs_Ven_AUC.svg']);
saveas(h5,[storage descr '_Art_vs_Ven_AUC.fig']);
            
figure(h6); supertitle('Arteriole vs Venule Peak Act: Intensity')
saveas(h6,[storage descr '_Art_vs_Ven_peak.svg']);
saveas(h6,[storage descr '_Art_vs_Ven_peak.fig']);

figure(h7); supertitle('Arteriole vs Venule Pixel Count: Intensity')
saveas(h7,[storage descr '_Art_vs_Ven_Ct.svg']);
saveas(h7,[storage descr '_Art_vs_Ven_Ct.fig']);

figure(h8); title('RM1 Full Timeseries')
saveas(h8,[storage descr '_RM1_Full_Timeseries_fUSI.svg']);
saveas(h8,[storage descr '_RM1_Full_Timeseries_fUSI.fig']);



%%
% figure;
% subplot(2,2,1); hold on
% stdshade(TMP{1,1}{5}',0.25,redC(1,:))
% stdshade(TMP{2,1}{5}',0.25,redC(2,:))
% stdshade(TMP{3,1}{5}',0.25,redC(3,:))
% subplot(2,2,2); hold on
% stdshade(TMP{1,1}{14}',0.25,redC(1,:))
% stdshade(TMP{2,1}{14}',0.25,redC(2,:))
% stdshade(TMP{3,1}{14}',0.25,redC(3,:))
% subplot(2,2,3); hold on
% stdshade(TMP{1,2}{5}',0.25,blueC(1,:))
% stdshade(TMP{2,2}{5}',0.25,blueC(2,:))
% stdshade(TMP{3,2}{5}',0.25,blueC(3,:))
% subplot(2,2,4); hold on
% stdshade(TMP{1,2}{14}',0.25,blueC(1,:))
% stdshade(TMP{2,2}{14}',0.25,blueC(2,:))
% stdshade(TMP{3,2}{14}',0.25,blueC(3,:))
% supertitle('M1 vascular')