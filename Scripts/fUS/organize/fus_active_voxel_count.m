function fus_active_voxel_count(storage,base_fold,slash,param)

% Extract number of significantly active voxels per ROI
roifolder = [storage(1:end-4) 'segment_ofUS' slash 'fus_20191205' slash];
roi = dir(roifolder);
roi(contains({roi.name},'.fig')) = [];
roi(~contains({roi.name},'R')) = [];

stim = {'0_1' '0_5' '1_0'};
descr = param.descr;

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
            cdata = imresize(ROI.cdata,.5);
        
            for i_mouse = 1:size(base_fold,1)
                
                pts{i_roi,i_mouse} = find(cdata); % indices of image in roi
                t_stat_mouse = t_stat(:,:,i_mouse);
                t_stat_roi{i_roi,i_mouse} = mean(t_stat_mouse(pts{i_roi,i_mouse})); % Average t-stat in roi
                
                p2tail_mouse = p2tail(:,:,i_mouse);
                p2tail_roi{i_roi,i_mouse} = p2tail_mouse(pts{i_roi,i_mouse});
                sig_pts_roi{i_roi,i_mouse} = sum(p2tail_roi{i_roi,i_mouse} < param.Pthresh); % number of sig active pixels
                sig_prop_roi{i_roi,i_mouse} = sig_pts_roi{i_roi,i_mouse}/size(pts{i_roi,i_mouse},1); % proportion of sig active pixels
                
                % FLOW ANALYSIS
                flow_file = [storage base_fold{i_mouse}(1:8) slash base_fold{i_mouse}(10:end) slash 'Flowc_ds.mat'];
                if exist(flow_file,'file')
                    
                    load(flow_file)
                    flow_roi{i_roi,i_mouse} = Fcds.signed(pts{i_roi,i_mouse}); 
                    veloc_roi{i_roi,i_mouse} = Fcds.veloc(pts{i_roi,i_mouse});
%                     if i_stim == 1 && i_roi == 1
%                         figure; imagesc(Fcds.veloc)
%                         map = Fcds.signed;
%                         map(map<=0) = 0; map(map>0)=1;
%                         figure; imagesc(map)
%                         cmap = ([1 1 1;60/255 72/255 142/255]); colormap(cmap)
%                         map = Fc.signed;
%                         map(map>0) = 0; map(map<0)=1;
%                         figure; imagesc(map)
%                         cmap = ([1 1 1;205/255 32/255 41/255]); colormap(cmap)
%                     end
                    
                    % Pos vs Neg flow: Total
                    flow_roi_pos_pts{i_roi,i_mouse} = find(flow_roi{i_roi,i_mouse} > 0); 
                    flow_roi_neg_pts{i_roi,i_mouse} = find(flow_roi{i_roi,i_mouse} < 0);
                    
                    % Pos vs Neg flow: Active
                    flow_roi_pos_sig_pts{i_roi,i_mouse} = find(flow_roi{i_roi,i_mouse} > 1 & p2tail_roi{i_roi,i_mouse} < param.Pthresh); 
                    flow_roi_neg_sig_pts{i_roi,i_mouse} = find(flow_roi{i_roi,i_mouse} < -1 & p2tail_roi{i_roi,i_mouse} < param.Pthresh);
                    
                    % Slow vs fast: Total
                    flow_roi_fast_pts{i_roi,i_mouse} = find(abs(veloc_roi{i_roi,i_mouse}) > 2.5);
                    flow_roi_slow_pts{i_roi,i_mouse} = find(abs(veloc_roi{i_roi,i_mouse}) < 2.5);
                    
                    % Slow vs fast: Active
                    flow_roi_fast_sig_pts{i_roi,i_mouse} = find(abs(veloc_roi{i_roi,i_mouse}) > 2.5 & p2tail_roi{i_roi,i_mouse} < param.Pthresh); 
                    flow_roi_slow_sig_pts{i_roi,i_mouse} = find(abs(veloc_roi{i_roi,i_mouse}) < 2.5 & p2tail_roi{i_roi,i_mouse} < param.Pthresh);

                end
            end
        end
        
        fprintf('\nSaving: %s\n', save_file);
        save(save_file,'pts','t_stat_roi','p2tail_roi','sig_pts_roi',...
            'sig_prop_roi','flow_roi','flow_roi_pos_pts','flow_roi_neg_pts',...
            'flow_roi_pos_sig_pts','flow_roi_neg_sig_pts','flow_roi_fast_pts',...
            'flow_roi_slow_pts','flow_roi_fast_sig_pts','flow_roi_slow_sig_pts');
    end
    
end
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                