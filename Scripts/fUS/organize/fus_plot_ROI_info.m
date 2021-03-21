function fus_plot_ROI_info(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};

fusroi = param.fusroi;
descr = param.descr;

h1 = figure(1);
h2 = figure(2);

clear StatROIt StatROIprop
StatROIt(1,1:3) = {'subjid','intensity','pairInt'};
id = repmat(1:size(base_fold,1),3,1); StatROIt(2:1+size(base_fold,1)*3,1) = num2cell(id(:));
StatROIt(2:1+size(base_fold,1)*3,2) = num2cell(repmat([1;2;3],size(base_fold,1),1));
StatROIt(2:1+size(base_fold,1)*3,3) = num2cell(repmat([1;2;3],size(base_fold,1),1));
StatROIprop = StatROIt;

for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_ind_fixed_effects.mat'];
    count_file = [stim_storage stim{i_stim} '_active_voxel_count.mat'];
    
    if exist(stim_file,'file') && exist(count_file,'file')
        
        load(stim_file)
        load(count_file)
        for i_roi = 1:20
            
            % Proportion ROI activated
            clear prop_roi
            for i_mouse = 1:size(base_fold,1)
                prop_roi(i_mouse) = size(find(p2tail_roi{i_roi,i_mouse} < param.Pthresh),1)/...
                    size(p2tail_roi{i_roi,i_mouse},1);
            end
            prop_roi_mean{i_roi,i_stim} = nanmean(prop_roi);
            prop_roi_sem{i_roi,i_stim} = nanstd(prop_roi)/sum(size(~isnan(prop_roi),2));
            
            % Format proportion activated values for stat testing (R .txt file)
            StatROIprop{1,3+i_roi} = fusroi{i_roi};
            StatROIprop(i_stim+1:3:1+size(base_fold,1)*3,3+i_roi) = num2cell(prop_roi');
            
            % Average T-val ROI activated
            tval_roi_mean{i_roi,i_stim} = nanmean(vertcat(t_stat_roi{i_roi,:}));
            tval_roi_sem{i_roi,i_stim} = nanstd(vertcat(t_stat_roi{i_roi,:}))/...
                sqrt(sum(~isnan(vertcat(t_stat_roi{i_roi,:}))));
            
            % Format average t-val values for stat testing (R .txt file)
            StatROIt{1,3+i_roi} = fusroi{i_roi};
            StatROIt(i_stim+1:3:1+size(base_fold,1)*3,3+i_roi) = num2cell(vertcat(t_stat_roi{i_roi,:}));
            
        end
    end
end

if isunix
    Tprop = ['/home/bradley/Dropbox/fUSI_ROI_prop_' descr '.txt'];
    Tt = ['/home/bradley/Dropbox/fUSI_ROI_Tval_' descr '.txt'];
elseif ispc
    Tprop = ['C:\Users\Brad\Dropbox\fUSI_ROI_prop_' descr '.txt'];
    Tt = ['C:\Users\Brad\Dropbox\fUSI_ROI_Tval_' descr '.txt'];
end

% Write stat file for R - ROI proportion activated
StatROIprop(:,[9,23]) = []; % Remove S2 regions!
fid = fopen(Tprop,'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatROIprop{1,:});
for K = 2:size(StatROIprop,1)
    fprintf(fid, '%.0f %.0f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatROIprop{K,:});
end
fclose(fid)

% Write stat file for R - Average ROI T-val
StatROIt(:,[9,23]) = []; % Remove S2 regions!
fid = fopen(Tt,'w');
fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', StatROIt{1,:});
for K = 2:size(StatROIt,1)
    fprintf(fid, '%.0f %.0f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', StatROIt{K,:});
end
fclose(fid)

figure(h1); clf
ctrs = 1:20; data = cell2mat(prop_roi_mean);
hBar{1} = bar(ctrs, data); hold on
clear ctr ydt
for k1 = 1:size(prop_roi_mean,2)
    ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
    ydt(k1,:) = hBar{1}(k1).YData;
end
errorbar(ctr, ydt, zeros(3,20), cell2mat(prop_roi_sem)', zeros(3,20), zeros(3,20), '.k',...
    'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
hBar{2} = bar(ctrs, data); delete(hBar{1})
set(gca,'ylim',[0 1],'xtick',1:20,'xticklabel',fusroi); ylabel('Proportion')

figure(h2); clf
ctrs = 1:20; data = cell2mat(tval_roi_mean);
hBar{1} = bar(ctrs, data); hold on
clear ctr ydt
for k1 = 1:size(tval_roi_mean,2)
    ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
    ydt(k1,:) = hBar{1}(k1).YData;
end
errorbar(ctr, ydt, zeros(3,20), cell2mat(tval_roi_sem)', zeros(3,20), zeros(3,20), '.k',...
    'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
hBar{2} = bar(ctrs, data); delete(hBar{1})
set(gca,'ylim',[0 11],'xtick',1:20,'xticklabel',fusroi); ylabel('Average T-val')
hold off; supertitle('Average ROI T-value')

figure(h1); supertitle('Proportion ROI Activated')
saveas(h1,[storage descr '_Prop_ROI_Act.svg']);
saveas(h1,[storage descr '_Prop_ROI_Act.fig']);

figure(h2); supertitle('Average ROI T-value')
saveas(h2,[storage descr '_Ave_ROI_T_val.svg']);
saveas(h2,[storage descr '_Ave_ROI_T_val.fig']);
