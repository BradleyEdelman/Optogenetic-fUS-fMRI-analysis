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

SingleROI = 'RM1';
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
        
        AUC_pos_tot{i_stim} = tsnormaveAUC_roi_sigpos;
        AUC_neg_tot{i_stim} = tsnormaveAUC_roi_signeg;
        CT_pos_tot{i_stim} = flow_roi_pos_sig_pts;
        CT_neg_tot{i_stim} = flow_roi_neg_sig_pts;
        
        for i_roi = 1:20
            
            % Anatomical ROI timeseries (across intensities)
            figure(h1); subplot(5,4,i_roi); hold on
            stdshade(vertcat(tsnormave_roi{i_roi,:}),0.25,stimC(i_stim,:))
            y1=get(gca,'ylim');
            plot([11 11],y1,'k'); plot([19 19],y1,'k')
            title(fusroi{i_roi})
            
            % Full timeseries for local stimulation site
            if strcmp(fusroi{i_roi},SingleROI)
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
                set(gca,'ylim', [-40 100]);
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
    
    % Plot noise z-score
    NOISE_Z{i_stim} = cell2mat(tsnormaveZ_noise(end,:));
    ROI_Z{i_stim} = tsnormaveZ_roi;
    
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

figure(h8); title('RM1 Full Timeseries')
saveas(h8,[storage descr '_' SingleROI '_Full_Timeseries_fUSI.svg']);
saveas(h8,[storage descr '_' SingleROI 'RM1_Full_Timeseries_fUSI.fig']);

if isunix
    Tvasauc = ['/home/bradley/Dropbox/fUSI_Vascular_AUC_' descr '.txt'];
    Tvaspk = ['/home/bradley/Dropbox/fUSI_Vascular_Pk_' descr '.txt'];
    Tvasct = ['/home/bradley/Dropbox/fUSI_Vascular_Ct_' descr '.txt'];
elseif ispc
    Tvasauc = ['C:\Users\bedelman\Dropbox\fUSI_Vascular_AUC_' descr '.txt'];
    Tvaspk = ['C:\Users\bedelman\Dropbox\fUSI_Vascular_Pk_' descr '.txt'];
    Tvasct = ['C:\Users\bedelman\Dropbox\fUSI_Vascular_Ct_' descr '.txt'];
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

% same as above but boxplots with data points
h9 = figure(9); h10 = figure(10);
num_an = size(tsnorm_roi,2);
for i_roi = 1:20
    
    % Plot AUC for arterial and venous components
    figure(h9); subplot(5,4,i_roi); cla; hold on
    D = [cell2mat(AUC_pos_tot{1}(i_roi,:))',...
        cell2mat(AUC_pos_tot{2}(i_roi,:))',...
        cell2mat(AUC_pos_tot{3}(i_roi,:))',...
        zeros(num_an,1),...
        cell2mat(AUC_neg_tot{1}(i_roi,:))',...
        cell2mat(AUC_neg_tot{2}(i_roi,:))',...
        cell2mat(AUC_neg_tot{3}(i_roi,:))'];
    Dl = [repmat('1',num_an,1);...
        repmat('2',num_an,1);...
        repmat('3',num_an,1);...
        repmat('4',num_an,1);...
        repmat('5',num_an,1);...
        repmat('6',num_an,1);...
        repmat('7',num_an,1)];
    cla
    boxplot(D,Dl); hold on
    set(gca,'xticklabel',{'0.1' '0.5' '1.0' '0' '0.1' '0.5' '1.0'})
    Xr = ones(size(D)).*(1+(rand(size(D))-0.5)/5);
    Xr(:,2) = Xr(:,2) + 1; Xr(:,3) = Xr(:,3) + 2;
    Xr(:,4) = Xr(:,4) + 3; Xr(:,5) = Xr(:,5) + 4;
    Xr(:,6) = Xr(:,6) + 5; Xr(:,7) = Xr(:,7) + 6;
    scatter(Xr(1:num_an*4),D(1:num_an*4),5,'r','filled')
    scatter(Xr(num_an*4+1:end),D(num_an*4+1:end),5,'b','filled')
    title(fusroi{i_roi})
    
    % Plot AUC for arterial and venous components
    figure(h10); subplot(5,4,i_roi); cla; hold on
    cellsz1 = cell2mat(cellfun(@size,CT_pos_tot{1}(i_roi,:),'uni',false));
    cellsz2 = cell2mat(cellfun(@size,CT_pos_tot{2}(i_roi,:),'uni',false));
    cellsz3 = cell2mat(cellfun(@size,CT_pos_tot{2}(i_roi,:),'uni',false));
    cellsz4 = cell2mat(cellfun(@size,CT_neg_tot{1}(i_roi,:),'uni',false));
    cellsz5 = cell2mat(cellfun(@size,CT_neg_tot{2}(i_roi,:),'uni',false));
    cellsz6 = cell2mat(cellfun(@size,CT_neg_tot{3}(i_roi,:),'uni',false));
    
    D11 = [cellsz1(1:2:end)',cellsz2(1:2:end)',cellsz3(1:2:end)',...
        zeros(num_an,1),...
        cellsz4(1:2:end)',cellsz5(1:2:end)',cellsz6(1:2:end)'];
    Dl2 = [repmat('1',num_an,1); repmat('2',num_an,1);...
        repmat('3',num_an,1); repmat('4',num_an,1);...
        repmat('5',num_an,1); repmat('6',num_an,1);...
        repmat('7',num_an,1)];
    cla
    boxplot(D11,Dl2); hold on
    set(gca,'xticklabel',{'0.1' '0.5' '1.0' '0' '0.1' '0.5' '1.0'})
    Xr = ones(size(D11)).*(1+(rand(size(D11))-0.5)/5);
    Xr(:,2) = Xr(:,2) + 1; Xr(:,3) = Xr(:,3) + 2;
    Xr(:,4) = Xr(:,4) + 3; Xr(:,5) = Xr(:,5) + 4;
    Xr(:,6) = Xr(:,6) + 5; Xr(:,7) = Xr(:,7) + 6;
    scatter(Xr(1:num_an*4),D11(1:num_an*4),5,'r','filled')
    scatter(Xr(num_an*4+1:end),D11(num_an*4+1:end),5,'b','filled')
    title(fusroi{i_roi})
    
    % Plot peak time series voxel-wise for ROIs
    if contains(['RM1' 'LM1' 'RCPu' 'LCPu'] ,fusroi{i_roi})
        f1 = figure;
        subplot(1,2,1); boxplot(D,Dl);
        set(gca,'xticklabel',{'0.1P' '0.5' '1.0' '0' '0.1' '0.5' '1.0'})
        title('AUC')
        if contains(['RM1', 'LM1'], fusroi{i_roi})
            set(gca,'ylim',[0 3500]);
        else
            set(gca,'ylim',[0 1500]);
        end
        
        subplot(1,2,2); boxplot(D11,Dl2);
        set(gca,'xticklabel',{'0.1N' '0.5' '1.0' '0' '0.1' '0.5' '1.0'})
        title('Count')
        if contains(['RM1', 'LM1'], fusroi{i_roi})
            set(gca,'ylim',[0 150]);
        else
            set(gca,'ylim',[0 475]);
        end
        
        saveas(f1,[storage descr '_AUC_box_' fusroi{i_roi} '.svg']);
    end    
    
end
figure(h9); supertitle('Arteriole vs Venule AUC: Boxplots')
saveas(h9,[storage descr '_AUC_box.svg']);
saveas(h9,[storage descr '_AUC_box.fig']);
figure(h10); supertitle('Arteriole vs Venule Act Count: Boxplots')
saveas(h10,[storage descr '_CT_box.svg']);
saveas(h10,[storage descr '_CT_box.fig']);

%% ROI Z SCORES

roi_int = [];
roi_int = [roi_int find(contains(fusroi, 'RM1'))];
roi_int = [roi_int find(contains(fusroi, 'LM1'))];
roi_int = [roi_int find(contains(fusroi, 'RCPu'))];
roi_int = [roi_int find(contains(fusroi, 'LCPu'))];

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
    title(fusroi{roi_int(i_roi)})
end
saveas(f,[storage descr '_ROI_CIRCUIT_Z.fig']);
saveas(f,[storage descr '_ROI_CIRCUIT_Z.svg']);

%% NOISE Z SCORES
% PLOT NOISE ROI Z-SCORE INFO
h27 = figure;
subplot(2,1,1);
D = [NOISE_Z{1}' NOISE_Z{2}' NOISE_Z{3}'];
Dl = ['1'; '1'; '1'; '1'; '1'; '1'; '1';...
    '2'; '2'; '2'; '2'; '2'; '2'; '2';...
    '3'; '3'; '3'; '3'; '3'; '3'; '3'];
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