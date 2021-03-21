function fus_regression_analysis(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};

fusroi = param.fusroi;
descr = param.descr;


sig_vein_t = cell(1,3);
sig_art_t = cell(1,3);
for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_ind_fixed_effects.mat'];
    count_file = [stim_storage stim{i_stim} '_active_voxel_count.mat'];
    ind_glm_file = [stim_storage stim{i_stim} '_ind_fixed_effects.mat'];
    
    if exist(stim_file,'file') && exist(count_file,'file') && exist(ind_glm_file,'file')
        
        load(count_file)
        load(stim_file)
        load(ind_glm_file)
        tval(:,:,i_stim) = reshape(t_stat,size(t_stat,1) * size(t_stat,2),size(t_stat,3));
%         pval(:,:,i_stim) = reshape(p2tail,size(t_stat,1) * size(t_stat,2),size(t_stat,3));
        
        sig_vein_vox = cell(1,size(base_fold,1));
        sig_art_vox = cell(1,size(base_fold,1));
        for i_roi = 1:10
            roi_pts = pts{i_roi,1}; % extract absolute voxel numbers in each roi
            
            for i_mouse = 1:size(base_fold,1)
                sig_vein_vox{i_mouse} = [sig_vein_vox{i_mouse}; roi_pts(flow_roi_pos_sig_pts{i_roi,i_mouse})]; % extract sig active veins from each roi
                sig_art_vox{i_mouse} = [sig_art_vox{i_mouse}; roi_pts(flow_roi_neg_sig_pts{i_roi,i_mouse})];
            end
        end
        
        
        % Extract t values for each mouse for vein vs art
        for i_mouse = 1:size(base_fold,1)
            
            sig_vein_vox_mouse = vertcat(sig_vein_vox{:,i_mouse});
            sig_vein_t{i_stim} = [sig_vein_t{i_stim}; tval(sig_vein_vox_mouse,i_mouse,i_stim)];
            
            sig_art_vox_mouse = vertcat(sig_art_vox{:,i_mouse});
            sig_art_t{i_stim} = [sig_art_t{i_stim}; tval(sig_art_vox_mouse,i_mouse,i_stim)];
            
        end
        
    end
end

for i = 1:3
    Vm(i) = mean(sig_vein_t{i});
    Ve(i) = std(sig_vein_t{i})/sqrt(size(sig_vein_t{i},1));
    Am(i) = mean(sig_art_t{i});
    Ae(i) = std(sig_art_t{i})/sqrt(size(sig_art_t{i},1));
end

figure; hold on
errorbar(Vm,Ve);
errorbar(Am,Ae)



% % REGRESSION FOR All ACTIVE VOXEL
% c_roi = [2:5,7:14,16:19];
% cortex_vein_t = cell(1,3); cortex_art_t = cell(1,3);
% for i_mouse = 1:size(base_fold,1)
%     for i_roi = 1:20
%         
%     
%     
%     cortex_vein = vertcat(flow_roi_pos_sig_pts{c_roi,i_mouse});
%     cortex_art = vertcat(flow_roi_neg_sig_pts{c_roi,i_mouse});
%     
%     for i_stim = 1:3
%         cortex_vein_t{i_stim} = [cortex_vein_t{i_stim}; tval(cortex_vein,i_mouse,i_stim)];
%         cortex_art_t{i_stim} = [cortex_art_t{i_stim}; tval(cortex_art,i_mouse,i_stim)];
%     end
% end
% 
% T1 = vertcat(cortex_vein_t{:});
% T2 = [0.1*ones(size(cortex_vein_t{1},1),1);0.5*ones(size(cortex_vein_t{1},1),1);1.0*ones(size(cortex_vein_t{1},1),1)];
% T2(isnan(T1)) = [];
% T1(isnan(T1)) = [];
% [B,BINT,R,RINT,STATS] = regress(T1,[ones(size(T2,1),1),T2]);
% ytmp = reshape(T1,[],3);    
% xtmp = repmat([0.1 0.5 1.0],[size(ytmp,1) 1]);
% 
% figure; cla; hold on
% scatter(xtmp(:),ytmp(:),10,'k','filled','marker','d')
% best = [0.1 0.5 1.0]*B(1) + B(2);
% plot([0.1 0.5 1.0],best,'b')
% % text(-.25,14,['y = ' num2str(B_roi(i)) '*b +' num2str(B_roi_int(i))]);
% set(gca,'ylim',[0 20],'xlim',[0 1.1])
% 
% T1 = vertcat(cortex_art_t{:});
% T2 = [0.1*ones(size(cortex_art_t{1},1),1);0.5*ones(size(cortex_art_t{1},1),1);1.0*ones(size(cortex_art_t{1},1),1)];
% T2(isnan(T1)) = [];
% T1(isnan(T1)) = [];
% [B,BINT,R,RINT,STATS] = regress(T1,[ones(size(T2,1),1),T2]);
% ytmp = reshape(T1,[],3);    
% xtmp = repmat([0.1 0.5 1.0],[size(ytmp,1) 1]);
% 
% figure; cla; hold on
% scatter(xtmp(:),ytmp(:),10,'k','filled','marker','d')
% best = [0.1 0.5 1.0]*B(1) + B(2);
% plot([0.1 0.5 1.0],best,'b')
% % text(-.25,14,['y = ' num2str(B_roi(i)) '*b +' num2str(B_roi_int(i))]);
% set(gca,'ylim',[0 20],'xlim',[0 1.1])
% 
% for i_pix = 1:size(tval,1)
%     
%     T1 = [tval(i_pix,:,1) tval(i_pix,:,2) tval(i_pix,:,3)]';
%     T2 = [0.1*ones(size(base_fold,1),1);0.5*ones(size(base_fold,1),1);1.0*ones(size(base_fold,1),1)];
% 
%     T2(isnan(T1)) = [];
%     T1(isnan(T1)) = [];
% 
%     [B,BINT,R,RINT,STATS] = regress(T1,[ones(size(T2,1),1),T2]);
%     
%     B_pix_int(i_pix) = B(1);
%     B_pix(i_pix) = B(2);
%     Rsq_pix(i_pix) = STATS(1);
%     P_pix(i_pix) = STATS(3);
%     
% end

% SINGLE VOXEL VEINS VS ARTERIES
c_roi = [2:5,7:14,16:19];
cortex_vein = vertcat(flow_roi_pos_sig_pts{c_roi},1);
cortex_art = vertcat(flow_roi_neg_sig_pts{c_roi},1);
vein_p = P_pix(cortex_vein);
art_p = P_pix(cortex_art);
vein_act = find(vein_p < 0.05);
art_act = find(art_p < 0.05);

% Plot veins and arteries
B_cvein = B_pix(cortex_vein); B_cvein_sig = B_cvein(vein_act);
B_int_cvein = B_pix(cortex_vein); B_int_cvein_sig = B_int_cvein(vein_act);

B_cart = B_pix(cortex_art); B_cart_sig = B_cart(art_act);
B_int_cart = B_pix(cortex_art); B_int_cart_sig = B_int_cart(art_act);

h5 = figure(25); subplot(1,3,1); hold on
scatter(B_cvein_sig,B_int_cvein_sig,5,'b','filled');
scatter(B_cart_sig,B_int_cart_sig,5,'b','filled');
set(gca,'xlim',[-10 10],'ylim',[-5 20]); title('fusi')
subplot(1,3,2); histogram(B_int_pix_sig,'BinEdges',[-5:1:10],'Normalization','probability');
set(gca,'xlim',[-5 10],'ylim',[0 .5]); title('int')
subplot(1,3,3); histogram(B_pix_sig,'BinEdges',[-5:1:10],'Normalization','probability');
set(gca,'xlim',[-5 15],'ylim',[0 .25]); title('beta')


cortex_fast = vertcat(flow_roi_fast_sig_pts{c_roi},1);
cortex_slow = vertcat(flow_roi_slow_sig_pts{c_roi},1);
fast_p = P_pix(cortex_fast);
slow_p = P_pix(cortex_slow);
fast_act = find(fast_p < 0.05);
slow_act = find(slow_p < 0.05);


% SINGLE VOXEL ANALYSIS
cortex_pix = vertcat(pts{[2:5,7:14,16:19],1}); % all cortical voxels

cortex_reg_p = P_pix(cortex_pix); % cortical voxel regression p values
cortex_reg_act = find(cortex_reg_p < 0.05); % 
cortex_prop_reg_sig = size(cortex_reg_act,2)/size(cortex_reg_p,2)*100 % proportion of cortical voxels sig linear 

glm_p = reshape(pval_grp,8800,3,1); % GLM p values
cortex_glm_p = glm_p(cortex_pix,:); % cortical GLM p values
cortex_glm_act = find(cortex_glm_p(:,1) < 0.05 &...
    cortex_glm_p(:,2) < 0.05 & cortex_glm_p(:,3) < 0.05); % find activated cortical voxels
cortex_prop_glm_sig = size(cortex_glm_act,1)/size(cortex_reg_p,2)*100 % proportion of cortical voxels active at all stim int

h4 = figure(77); histogram(cortex_reg_p,40,'Normalization','probability');
set(gca,'ylim',[0 .65]); title('fusi')

% Plot all significant brain pixels
B_pix_sig = B_pix(cortex_pix); B_pix_sig = B_pix_sig(cortex_reg_act);
B_int_pix_sig = B_pix_int(cortex_pix); B_int_pix_sig = B_int_pix_sig(cortex_reg_act);
h5 = figure(78); subplot(1,3,1); scatter(B_int_pix_sig,B_pix_sig,5,'b','filled');
set(gca,'xlim',[-10 10],'ylim',[-5 20]); title('fusi')
subplot(1,3,2); histogram(B_int_pix_sig,'BinEdges',[-5:1:10],'Normalization','probability');
set(gca,'xlim',[-5 10],'ylim',[0 .5]); title('int')
subplot(1,3,3); histogram(B_pix_sig,'BinEdges',[-5:1:10],'Normalization','probability');
set(gca,'xlim',[-5 15],'ylim',[0 .25]); title('beta')

% Plot all pixels in local region (M1)
% h6 = figure(79);
% for i_roi = 1:20
%     subplot(5,4,i_roi)
%     B_int_tmp = B_pix_int(pts{i_roi,1});
%     B_tmp = B_pix(pts{i_roi,1});
%     P_tmp = P_pix(pts{i_roi,1});
%     scatter(B_int_tmp(P_tmp < 0.05),B_tmp(P_tmp < 0.05),5,'b','filled');
%     set(gca,'xlim',[-10 10],'ylim',[-5 20]);
%     title(param.fusroi{i_roi});
% end
% 
% h6 = figure(79); subplot(1,3,1); scatter(B_int_tmp(P_tmp < 0.05),B_tmp(P_tmp < 0.05),5,'b','filled');
% set(gca,'xlim',[-10 10],'ylim',[-5 20]); title('fusi')
% subplot(1,3,2); histogram(B_pix_int(pts{14,1}),'BinEdges',[-5:1:20],'Normalization','probability');
% set(gca,'xlim',[-5 10]); title('int')
% subplot(1,3,3); histogram(B_pix(pts{14,1}),'BinEdges',[-5:1:20],'Normalization','probability');
% set(gca,'xlim',[-5 15]); title('beta')


% ROI REGRESSION ANALYSIS USING ONLY SIG LINEAR VOXELS
h6 = figure(79);
for i = 1:20
    pts_roi = pts{i,1};
    
    T1_roi = [];
    count = 0;
    for j = 1:size(pts_roi,1)
        if P_pix(pts_roi(j)) < 0.05
            t_tmp = [tval(pts_roi(j),:,1) tval(pts_roi(j),:,2) tval(pts_roi(j),:,3)]';
            T1_roi = [T1_roi t_tmp];
            count = count + 1;
        end
    end
    T1_roi = nanmean(T1_roi,2);
    T2 = [0.1*ones(size(base_fold,1),1);0.5*ones(size(base_fold,1),1);1.0*ones(size(base_fold,1),1)];
    
    T2(isnan(T1_roi)) = [];
    T1_roi(isnan(T1_roi)) = [];
    if ~isempty(T1_roi)
        [B,BINT,R,RINT,STATS] = regress(T1_roi,[ones(size(T2,1),1),T2]);
        
        ytmp = reshape(T1_roi,[],3);    
        xtmp = repmat([0.1 0.5 1.0],[size(ytmp,1) 1]);
        
        subplot(5,4,i); cla; hold on
        scatter(xtmp(:),ytmp(:),10,'k','filled','marker','d')
        best = [0.1 0.5 1.0]*B_roi(i) + B_roi_int(i);
        
        plot([0.1 0.5 1.0],best,'b')
        text(-.25,14,['y = ' num2str(B_roi(i)) '*b +' num2str(B_roi_int(i))]);
        set(gca,'ylim',[0 20],'xlim',[0 1.1])
        title(fusroi{i})

        B_roi_int_sig(i) = B(1);
        B_roi_sig(i) = B(2);
        Rsq_roi_sig(i) = STATS(1);
        P_roi_sig(i) = STATS(3);
        count_roi_sig(i) = count;
        prop_roi_sig(i) = count/size(pts_roi,1);
        
    else
        
        B_roi_int_sig(i) = 0;
        B_roi_sig(i) = 0;
        Rsq_roi_sig(i) = 0;
        P_roi_sig(i) = 1;
        count_roi_sig(i) = 0;
        prop_roi_sig(i) = 0;
        
    end
end

B_pix_int = reshape(B_pix_int,110,80); B_pix = reshape(B_pix,110,80);
Rsq_pix = reshape(Rsq_pix,110,80); P_pix = reshape(P_pix,110,80);

% rescale dynamic range of anat for vis
maxR = max(Rsq_pix(:)); maxB = ceil(max(abs(B_pix(:)))); plotmin = -100; plotmax = 0; maxB = 15; maxB2 = 8;
CbarR = [gray(abs(plotmax - plotmin)*100); fireice(2*maxR*100)]; CAXR = [plotmin - maxR maxR];
CbarB = [gray(abs(plotmax - plotmin)*100); fireice(2*maxB*100)]; CAXB = [plotmin - maxB maxB];
CbarB2 = [gray(abs(plotmax - plotmin)*100); fireice(2*maxB2*100)]; CAXB2 = [plotmin - maxB2 maxB2];

atmp = [storage '20191119' slash '4346075_N_D2' slash 'Anat.mat'];
anat = load(atmp); anat = anat.ref;
BB22 = imresize(B_pix_int,size(anat));
BB2 = imresize(B_pix,size(anat));
Rsq2 = imresize(Rsq_pix,size(anat));
PP2 = imresize(P_pix,size(anat));
anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin; 
anatRR = anat - maxR; anatBB = anat - maxB; anatBB2 = anat - maxB2;
% anatBB(PP2 < .05) = BB2(PP2 < .05); anatBB(mask(:,:) == 0) = CAXB(1);
% anatRR(PP2 < .05) = RR2(PP2 < .05); anatRR(mask(:,:) == 0) = CAXR(1);

h3 = figure(3); clf
subplot(1,4,1); imagesc(Rsq2); title('Rsq pixel');
set(gca,'xtick',[],'ytick',[]); caxis([0 1]); colorbar
subplot(1,4,2); imagesc(PP2); title('P-val pixel');
set(gca,'xtick',[],'ytick',[]); caxis([0 0.05]); colorbar
subplot(1,4,3); imagesc(BB2); title('Beta pix');
set(gca,'xtick',[],'ytick',[]); caxis([0 maxB]); colorbar
subplot(1,4,4); imagesc(BB22); title('Int pix');
set(gca,'xtick',[],'ytick',[]); caxis([0 maxB2]); colorbar
colormap jet
supertitle('Intensity Regression: Pixel fUSI')

saveas(h1,[storage descr '_Intensity_Regression_ROI_Scatter_fUSI.svg']);
saveas(h1,[storage descr '_Intensity_Regression_ROI_Scatter_fUSI.fig']);

saveas(h2,[storage descr '_Intensity_Regression_ROI_fUSI.svg']);
saveas(h2,[storage descr '_Intensity_Regression_ROI_fUSI.fig']);

saveas(h3,[storage descr '_Intensity_Regression_Pix_fUSI.svg']);
saveas(h3,[storage descr '_Intensity_Regression_Pix_fUSI.fig']);

saveas(h4,[storage descr '_Intensity_Regression_Pix_Pval_hist_fUSI.svg']);
saveas(h4,[storage descr '_Intensity_Regression_Pix_Pval_hist_fUSI.fig']);

saveas(h5,[storage descr '_Intensity_Regression_Pix_beta_vs_int_fUSI.svg']);
saveas(h5,[storage descr '_Intensity_Regression_Pix_beta_vs_int_fUSI.fig']);

saveas(h6,[storage descr '_Intensity_Regression_ROI_sig_pix_fUSI.svg']);
saveas(h6,[storage descr '_Intensity_Regression_ROI_sig_pix_fUSI.fig']);

save_file = [storage descr '_regression.mat'];
fprintf('\nSaving: %s\n', save_file);
save(save_file,'B_roi_int','B_roi','Rsq_roi','P_roi','B_pix_int','B_pix','Rsq_pix','P_pix',...
    'cortex_prop_reg_sig','cortex_prop_glm_sig','B_pix_sig','B_int_pix_sig','B_roi_int_sig',...
    'B_roi_sig','Rsq_roi_sig','P_roi_sig','count_roi_sig','prop_roi_sig')