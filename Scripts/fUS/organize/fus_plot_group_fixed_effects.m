function fus_plot_group_fixed_effects(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
descr = param.descr;

% Load group data
for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_grp_data.mat'];
    stim_file = [stim_storage stim{i_stim} '_grp_fixed_effects.mat'];
    
    if exist(stim_file,'file')
        
        load(stim_file)
        
        % Plot Design Matrix
        h1 = figure(1); subplot(2,1,1);
        b=bar(.5:1:size(grp_c,1),grp_c,'k','barwidth',.5); grid minor;
        set(gca,'position',[.13 .80 .775 .15],'ylim',[-2 2],'xlim',[0 size(grp_c,1)],...
            'xtick',0:1:12,'xticklabel',[],'yticklabel',[])
        subplot(2,1,2); image(grp_D*64);
        set(gca,'position',[.13 .1 .775 .65],'xticklabel',[],'yticklabel',[1 size(grp_D,1)],'ytick',[1 size(grp_D,1)])
        ylabel('Acq Number');
        colormap gray
        supertitle(['Design Matrix: ' stim{i_stim}])
        
        % Create colorbar for activation maps
        plotmin = -100; plotmax = 0; maxt = param.maxt; Pthresh = param.Pthresh;
        Cbar = [gray(abs(plotmax - plotmin)*100); fireice(2*maxt*100)]; CAX = [plotmin - maxt maxt];
        Cbarp = [jet(round(1.1*100))];
        
        % Load anatomical image
        atmp = [storage '20191119' slash '4346075_N_D2' slash 'Anat.mat'];
        anii = load(atmp); anii = anii.ref;
        
        mask = ones(size(anii));
        % Resize activation to anatomical sampling
        grp_c_p2tail_tmp = imresize(grp_c_p2tail,[size(anii,1) size(anii,2)],'bilinear');
        grp_c_p2tail_tmp(mask == 0) = 1.1;
        grp_t_stat_tmp = imresize(grp_t_stat,[size(anii,1) size(anii,2)],'bilinear');
        grp_t_stat_tmp(mask == 0) = 0;

        % rescale dynamic range of anatomical for visuslization
        anat = anii;
        anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
        anat = anat - maxt;
        grp_anat = anat; grp_anat(grp_c_p2tail_tmp < Pthresh) = grp_t_stat_tmp(grp_c_p2tail_tmp < Pthresh); grp_anat(mask == 0) = 0;

        h2 = figure(2);
        subplot(size(stim,2),2,1+2*(i_stim-1))
        imagesc(grp_c_p2tail_tmp); set(gca,'xticklabel',[],'yticklabel',[]); colormap(gca,Cbarp); caxis([0 1.1]); title(['p-val 2tail - ' stim{i_stim} 'mW'])
        subplot(size(stim,2),2,2+2*(i_stim-1))
        imagesc(grp_anat); set(gca,'xticklabel',[],'yticklabel',[]); colormap(gca,Cbar); caxis(CAX); title(['t-val 2tail - ' stim{i_stim} 'mW'])

    end
end

saveas(h1,[stim_storage 'Grp_design_matrix.svg']);
saveas(h1,[stim_storage 'Grp_design_matrix.fig']);

saveas(h2,[stim_storage 'Grp_activation_maps.svg']);
saveas(h2,[stim_storage 'Grp_activation_maps.fig']);

