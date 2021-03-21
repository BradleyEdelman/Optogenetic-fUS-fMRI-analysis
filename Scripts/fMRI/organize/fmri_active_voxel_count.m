function fmri_active_voxel_count(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
descr = param.descr;

% Extract number of significantly active voxels per ROI
switch descr
    case 'pre' 
        roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_pre_20191205' slash];
    case 'single' 
        roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_pre_20191205' slash];
    case 'ephys'
        roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_pre_20191205' slash];
    case 'post' 
        roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_post_20200212' slash];
    case 'ctrl' 
        roifolder = [storage(1:end-5) 'segment_ofUS' slash 'fmri_ctrl_20200212' slash];
end
roi = dir(roifolder);
roi(contains({roi.name},'.fig')) = [];
roi(~contains({roi.name},'R')) = [];

for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_ind_fixed_effects.mat'];
    save_file = [stim_storage stim{i_stim} '_active_voxel_count.mat'];
    
    if exist(stim_file,'file')
        
        load(stim_file)
    
        for i_roi = 1:size(roi,1)
        
            roifile = [roifolder roi(i_roi).name];
            load(roifile)
            cdata = imresize(ROI.cdata,.2);
        
            for i_mouse = 1:size(base_fold,1)
                
                pts{i_roi,i_mouse} = find(cdata); % indices of image in roi
                t_stat_mouse = t_stat(:,:,i_mouse);
                t_stat_roi{i_roi,i_mouse} = mean(t_stat_mouse(pts{i_roi,i_mouse})); % Average t-stat in roi
                
                p2tail_mouse = p2tail(:,:,i_mouse);
                p2tail_roi{i_roi,i_mouse} = p2tail_mouse(pts{i_roi,i_mouse});
                sig_pts_roi{i_roi,i_mouse} = sum(p2tail_roi{i_roi,i_mouse} < param.Pthresh); % number of sig active pixels
                sig_prop_roi{i_roi,i_mouse} = sig_pts_roi{i_roi,i_mouse}/size(pts{i_roi,i_mouse},1); % proportion of sig active pixels
               
            end
        end
        
        fprintf('\nSaving: %s\n', save_file);
        save(save_file,'pts','t_stat_roi','p2tail_roi','sig_pts_roi','sig_prop_roi');
    end
    
end
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                