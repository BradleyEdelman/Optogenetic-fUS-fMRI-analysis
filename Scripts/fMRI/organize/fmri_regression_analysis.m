function fmri_regression_analysis(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};

fmriroi = param.fmriroi;
descr = param.descr;
sliceidx = param.sliceidx;

for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_ind_fixed_effects.mat'];
    count_file = [stim_storage stim{i_stim} '_active_voxel_count.mat'];
    grp_glm_file = [stim_storage stim{i_stim} '_grp_fixed_effects.mat'];
    
    if exist(stim_file,'file') && exist(count_file,'file')
        
        load(count_file)
        load(stim_file)
        load(grp_glm_file)
%         peak(:,:,i_stim) = cell2mat(tsnormavepeak);
%         peak_roi(:,:,i_stim) = cell2mat(tsnormavepeak_roi);
        tval(:,:,i_stim) = reshape(t_stat,size(t_stat,1) * size(t_stat,2),size(t_stat,3));
        pval_grp(:,:,i_stim) = reshape(grp_p2tail,size(t_stat,1) * size(t_stat,2),1);
        tval_roi(:,:,i_stim) = cell2mat(t_stat_roi);
        
    end
end


h1 = figure(1);
for i_roi = 1:20
    
    T1 = [tval_roi(i_roi,:,1) tval_roi(i_roi,:,2) tval_roi(i_roi,:,3)]';
    T2 = [0.1*ones(size(base_fold,1),1);0.5*ones(size(base_fold,1),1);1.0*ones(size(base_fold,1),1)];

    T2(isnan(T1)) = [];
    T1(isnan(T1)) = [];

    [B,BINT,R,RINT,STATS] = regress(T1,[ones(size(T2,1),1),T2]);
    
    B_roi_int(i_roi) = B(1);
    B_roi(i_roi) = B(2);
    Rsq_roi(i_roi) = STATS(1);
    P_roi(i_roi) = STATS(3);
    MINACT(i_roi) = B_roi_int(i_roi)+B_roi(i_roi)*0.1;
    
    subplot(5,4,i_roi); cla; hold on
    set(gca,'xlim',[-.5 1.5],'ylim',[-2.5 15]);
%     scatter(T2,T1); hold on
%     plot(reshape(T2,[],3)',reshape(T1,[],3)')
%     scatter([0.1 0.5 1.0], mean(reshape(T1,[],3),1),50,'k','filled','marker','sq')
%     plot([0.1 0.5 1.0],mean(reshape(T1,[],3),1),'k','linewidth',2)
    xtmp = repmat([0.1 0.5 1.0],[7 1]);
    ytmp = reshape(T1,[],3);
    scatter(xtmp(:),ytmp(:),10,'k','filled','marker','d')
    best = [0.1 0.5 1.0]*B_roi(i_roi) + B_roi_int(i_roi);
    plot([0.1 0.5 1.0],best,'b')
    text(-.25,14,['y = ' num2str(B_roi(i_roi)) '*b +' num2str(B_roi_int(i_roi))]);
    title(fmriroi{i_roi})
    
end

h2 = figure(2); clf
subplot(1,5,1); imagesc(Rsq_roi'); title('Rsq ROI');
set(gca,'ytick',1:20,'yticklabel',param.fmriroi,'xtick',[]); caxis([0 1]); colorbar
subplot(1,5,2); imagesc(P_roi'); title('P-val ROI');
set(gca,'ytick',1:20,'yticklabel',param.fmriroi,'xtick',[]); caxis([0 1]); colorbar
subplot(1,5,3); imagesc(B_roi'); title('Beta ROI');
set(gca,'ytick',1:20,'yticklabel',param.fmriroi,'xtick',[]); caxis([0 max(B_roi)]); colorbar
subplot(1,5,4); imagesc(B_roi_int'); title('Int ROI');
set(gca,'ytick',1:20,'yticklabel',param.fmriroi,'xtick',[]); caxis([0 max(B_roi_int)]); colorbar
subplot(1,5,5); imagesc(MINACT'); title('Model Predict 0.1mW');
set(gca,'ytick',1:20,'yticklabel',param.fmriroi,'xtick',[]); caxis([0 max(MINACT)]); colorbar
colormap jet
supertitle('Intensity Regression: ROI')

% REGRESSION FOR EACH VOXEL
for i_pix = 1:size(tval,1)
    
    T1 = [tval(i_pix,:,1) tval(i_pix,:,2) tval(i_pix,:,3)]';
    T2 = [0.1*ones(size(base_fold,1),1);0.5*ones(size(base_fold,1),1);1.0*ones(size(base_fold,1),1)];

    T2(isnan(T1)) = [];
    T1(isnan(T1)) = [];

    [B,BINT,R,RINT,STATS] = regress(T1,[ones(size(T2,1),1),T2]);
    
    B_pix_int(i_pix) = B(1);
    B_pix(i_pix) = B(2);
    Rsq_pix(i_pix) = STATS(1);
    P_pix(i_pix) = STATS(3);
    MINACT_pix(i_pix) = B_pix_int(i_pix)+B_pix(i_pix)*0.1;
    
end


% SINGLE VOXEL ANALYSIS
cortex_pix = vertcat(pts{[2:5,7:14,16:19],1}); % all cortical voxels

cortex_reg_p = P_pix(cortex_pix); % cortical voxel regression p values
cortex_reg_act = find(cortex_reg_p < 0.05); % 
cortex_prop_reg_sig = size(cortex_reg_act,2)/size(cortex_reg_p,2)*100 % proportion of cortical voxels sig linear 

glm_p = reshape(pval_grp,2800,3,1); % GLM p values
cortex_glm_p = glm_p(cortex_pix,:); % cortical GLM p values
cortex_glm_act = find(cortex_glm_p(:,1) < 0.05 &...
    cortex_glm_p(:,2) < 0.05 & cortex_glm_p(:,3) < 0.05); % find activated cortical voxels
cortex_prop_glm_sig = size(cortex_glm_act,1)/size(cortex_reg_p,2)*100 % proportion of cortical voxels active at all stim int

h4 = figure(77); histogram(cortex_reg_p,40,'Normalization','probability');
set(gca,'ylim',[0 .65]); title('fmri')

% Plot all significant brain pixels
B_pix_sig = B_pix(cortex_pix); B_pix_sig = B_pix_sig(cortex_reg_act);
B_int_pix_sig = B_pix_int(cortex_pix); B_int_pix_sig = B_int_pix_sig(cortex_reg_act);
h5 = figure(78); subplot(1,3,1); scatter(B_int_pix_sig,B_pix_sig,5,'b','filled');
set(gca,'xlim',[-10 10],'ylim',[-5 20]); title('fmri')
subplot(1,3,2); histogram(B_int_pix_sig,'BinEdges',[-5:1:10],'Normalization','probability');
set(gca,'xlim',[-5 10],'ylim',[0 .5]); title('int')
subplot(1,3,3); histogram(B_pix_sig,'BinEdges',[-5:1:10],'Normalization','probability');
set(gca,'xlim',[-5 15],'ylim',[0 .25]); title('beta')
supertitle('beta v intercept')

% Plot all significant brain pixels (PREDICT 0.1 vs SLOPE)
B_pix_sig = B_pix(cortex_pix); B_pix_sig = B_pix_sig(cortex_reg_act);
MINACT_pix_sig = MINACT_pix(cortex_pix); MINACT_pix_sig = MINACT_pix_sig(cortex_reg_act);
h6 = figure(89); subplot(1,3,1); scatter(MINACT_pix_sig,B_pix_sig,5,'b','filled');
set(gca,'xlim',[-10 10],'ylim',[-5 20]); title('fusi')
subplot(1,3,2); histogram(MINACT_pix_sig,'BinEdges',[-5:1:10],'Normalization','probability');
set(gca,'xlim',[-5 10],'ylim',[0 .5]); title('int')
subplot(1,3,3); histogram(B_pix_sig,'BinEdges',[-5:1:10],'Normalization','probability');
set(gca,'xlim',[-5 15],'ylim',[0 .25]); title('beta')
supertitle('beta vs model predict 0.1mW')

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
        best = [0.1 0.5 1.0]*B(2) + B(1);
        
        plot([0.1 0.5 1.0],best,'b')
        text(-.25,14,['y = ' num2str(B(2)) '*b +' num2str(B(1))]);
        set(gca,'ylim',[0 20],'xlim',[0 1.1])
        title(fmriroi{i})

        B_roi_int_sig(i) = B(1);
        B_roi_sig(i) = B(2);
        Rsq_roi_sig(i) = STATS(1);
        P_roi_sig(i) = STATS(3);
        MINACT_roi_sig(i) = B_roi_int_sig(i) + B_roi_sig(i)*0.1;
        count_roi_sig(i) = count;
        prop_roi_sig(i) = count/size(pts_roi,1);
        
    else
        
        B_roi_int_sig(i) = 0;
        B_roi_sig(i) = 0;
        Rsq_roi_sig(i) = 0;
        P_roi_sig(i) = 1;
        MINACT_roi_sig(i) = 0;
        count_roi_sig(i) = 0;
        prop_roi_sig(i) = 0;
        
    end
end

B_pix_int = reshape(B_pix_int,35,80); B_pix = reshape(B_pix,35,80);
Rsq_pix = reshape(Rsq_pix,35,80); P_pix = reshape(P_pix,35,80);
MINACT_pix = reshape(MINACT_pix,35,80);

% rescale dynamic range of anat for vis
maxR = max(Rsq_pix(:)); maxB = ceil(max(abs(B_pix(:)))); plotmin = -100; plotmax = 0; maxB = 15; maxB2 =8; maxB3 = 10;
CbarR = [gray(abs(plotmax - plotmin)*100); fireice(2*maxR*100)]; CAXR = [plotmin - maxR maxR];
CbarB = [gray(abs(plotmax - plotmin)*100); fireice(2*maxB*100)]; CAXB = [plotmin - maxB maxB];
CbarB2 = [gray(abs(plotmax - plotmin)*100); fireice(2*maxB2*100)]; CAXB2 = [plotmin - maxB2 maxB2];
CbarB3 = [gray(abs(plotmax - plotmin)*100); fireice(2*maxB3*100)]; CAXB3 = [plotmin - maxB3 maxB3];

anat = load_nii(param.atmp); anat = anat.img(:,:,sliceidx);
mask = load_nii(param.mtmp); mask = mask.img(:,:,sliceidx);
BB22 = imresize(B_pix_int,size(anat));
BB2 = imresize(B_pix,size(anat));
Rsq2 = imresize(Rsq_pix,size(anat));
PP2 = imresize(P_pix,size(anat));
MA2 = imresize(MINACT_pix, size(anat));

anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin; 
anatRR = anat - maxR; anatBB = anat - maxB; anatBB2 = anat - maxB2; anatPP = anat;
BB22(mask == 0) = CAXB2(1); BB2(mask == 0) = CAXB(1); Rsq2(mask == 0) = CAXR(1); PP2(mask == 0) = 1; MA2(mask == 0) = 0;

h3 = figure(3); clf
subplot(1,5,1); imagesc(Rsq2); title('Rsq pixel');
set(gca,'xtick',[],'ytick',[]); caxis([0 1]); colorbar
subplot(1,5,2); imagesc(PP2); title('P-val pixel');
set(gca,'xtick',[],'ytick',[]); caxis([0 .05]); colorbar
subplot(1,5,3); imagesc(BB2); title('Beta pix');
set(gca,'xtick',[],'ytick',[]); caxis([0 maxB]); colorbar
subplot(1,5,4); imagesc(BB22);title('Int pix');
set(gca,'xtick',[],'ytick',[]); caxis([0 maxB2]); colorbar
subplot(1,5,5); imagesc(MA2); title('model 0.1mW');
set(gca,'xtick',[],'ytick',[]); caxis([0 maxB3]); colorbar
colormap jet
supertitle('Intensity Regression: Pixel')

saveas(h1,[storage descr '_Intensity_Regression_ROI_Scatter_fMRI.svg']);
saveas(h1,[storage descr '_Intensity_Regression_ROI_Scatter_fMRI.fig']);

saveas(h2,[storage descr '_Intensity_Regression_ROI_fMRI.svg']);
saveas(h2,[storage descr '_Intensity_Regression_ROI_fMRI.fig']);

saveas(h3,[storage descr '_Intensity_Regression_Pix_fMRI.svg']);
saveas(h3,[storage descr '_Intensity_Regression_Pix_fMRI.fig']);

saveas(h4,[storage descr '_Intensity_Regression_Pix_Pval_hist_fUSI.svg']);
saveas(h4,[storage descr '_Intensity_Regression_Pix_Pval_hist_fUSI.fig']);

saveas(h5,[storage descr '_Intensity_Regression_Pix_beta_vs_int_fUSI.svg']);
saveas(h5,[storage descr '_Intensity_Regression_Pix_beta_vs_int_fUSI.fig']);

saveas(h6,[storage descr '_Intensity_Regression_ROI_sig_pix_fUSI.svg']);
saveas(h6,[storage descr '_Intensity_Regression_ROI_sig_pix_fUSI.fig']);

save_file = [storage descr '_regression.mat'];
fprintf('\nSaving: %s\n', save_file);
save(save_file,'MINACT', 'MINACT_pix', 'MINACT_pix_sig', 'MINACT_roi_sig', 'Rsq_roi','P_roi','B_pix_int','B_pix','Rsq_pix','P_pix',...
    'cortex_prop_reg_sig','cortex_prop_glm_sig','B_pix_sig','B_int_pix_sig','B_roi_int_sig',...
    'B_roi_sig','Rsq_roi_sig','P_roi_sig','count_roi_sig','prop_roi_sig')




