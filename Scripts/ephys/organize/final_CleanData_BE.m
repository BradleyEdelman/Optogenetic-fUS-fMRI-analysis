function final_CleanData_BE(storage,base_fold,slash,type)
% Define bad channels/trials/units to remove from the data

%%

stim = {'0_1' '0_5' '1_0'};

badChannels = []; badTrials = []; badUnits = [];
for i_mouse = 1:size(base_fold,1)
    
    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '_adj.mat'];
        clipidx_file = [data_fold stim{i_stim} '_clipidx.mat'];
        save_file = [data_fold stim{i_stim} '_clean.mat'];
    
        if exist(data_file,'file') && exist(clipidx_file,'file')
            % Load the data
            load(data_file);
            load(clipidx_file);
            
            if strcmp(type,'Thy1')
                if i_mouse == 4 && i_stim == 1
                    clipidx = [clipidx; 7*ones(32,1) (1:32)'; 8*ones(32,1) (1:32)'; 9*ones(32,1) (1:32)'];
                elseif i_mouse == 4 && i_stim == 2
                    clipidx = [clipidx; 1*ones(32,1) (1:32)'; 2*ones(32,1) (1:32)'; 6*ones(32,1) (1:32)'];
                elseif i_mouse == 1
                    clipidx = [clipidx; [1:10]' 11*ones(10,1); [1:10]' 12*ones(10,1); [1:10]' 13*ones(10,1);...
                        [1:10]' 14*ones(10,1); [1:10]' 15*ones(10,1); [1:10]' 16*ones(10,1); [1:10]' 28*ones(10,1);...
                        [1:10]' 29*ones(10,1); [1:10]' 30*ones(10,1); [1:10]' 31*ones(10,1); [1:10]' 32*ones(10,1);];
                elseif i_mouse == 1 && i_stim == 1
                    clipidx = [clipidx; 3*ones(32,1) (1:32)'];
                end
            elseif strcmp(type,'ctrl')
            end
            
            % Remove flagged trials and channels from LFP
            for i_clip = 1:size(clipidx,1)
                LFP_mV_adj{clipidx(i_clip,1),clipidx(i_clip,2)} = [];
                LFP_TS_adj{clipidx(i_clip,1),clipidx(i_clip,2)} = [];
                Laser_TS_adj{clipidx(i_clip,1),clipidx(i_clip,2)} = [];
                Spike_TS_adj{clipidx(i_clip,1),clipidx(i_clip,2)} = [];
            end
    
            % Save the data to a new file with the filters saved as well
            fprintf('Saving: %s\n\n',save_file);
            save(save_file,'LFP_TS_adj','LFP_mV_adj','Spike_TS_adj','Laser_TS_adj','badChannels','badTrials','badUnits');
        end
    end
end