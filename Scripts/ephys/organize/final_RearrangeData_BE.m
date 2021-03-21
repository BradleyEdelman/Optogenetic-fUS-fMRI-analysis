function final_RearrangeData_BE(storage,base_fold,slash)
% Load and crop the data into 20 trials of 20" on, 20" on, 20" off
% Also compare the identified units across different recordings and
% frequencies, make sure you include them all

%%
stim = {'0_1' '0_5' '1_0'};

for i_mouse = 1:size(base_fold,1)
    
    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '.mat'];
        save_file = [data_fold stim{i_stim} '_adj.mat'];
        
        if exist(data_file,'file')
            load(data_file)
            
            TS_delta = Laser_TS(1) - 24;
            if TS_delta < 0
                Diff = find(diff(Laser_TS) > 1);
                Diff = Diff(1) + 1;
                TS_delta = Laser_TS(Diff) - 24;
            end
            
            % Create a cell array to hold the shifted Laser_TS values for each trial
            temp_TS = Laser_TS - TS_delta;
            trialcnt = floor(max(temp_TS/60));
            for k = 1:trialcnt
                trialidx{k} = find(temp_TS >= 0 + (k-1)*60 & temp_TS < 60 + (k-1)*60);
                Laser_TS_adj{k} = temp_TS(trialidx{k}) - 60*(k-1);
            end                 
                        
            % Create a cell array to hold the shifted LFP_TS values for each trial
            temp_TS = LFP_TS - TS_delta;
            clear trialidx
            for k = 1:trialcnt
                trialidx{k} = find(temp_TS >= 0 + (k-1)*60 & temp_TS < 60 + (k-1)*60);
                LFP_TS_adj{k} = temp_TS(trialidx{k}) - 60*(k-1);
            end     
            
            % Create a cell array for LFP_mV aligned to the laser with all trials
            for i_channel = 1:size(LFP_mV,2)
                temp_mV = LFP_mV{i_channel};
                for i_trial = 1:trialcnt
                    LFP_mV_adj{i_trial,i_channel} = temp_mV(trialidx{i_trial});
                end
            end
                        
            % Create a Trial x Channel x Unit cell array for Spike_TS aligned to the laser
            for i_channel = 1:size(Spike_TS,2)
                for i_unit = 1:size(Spike_TS,3)
                    temp_TS = Spike_TS{1,i_channel,i_unit} - TS_delta;
                    
                    clear trialidx
                    for i_trial = 1:trialcnt
                        trialidx{i_trial} = find(temp_TS >= 0 + (i_trial-1)*60 & temp_TS < 60 + (i_trial-1)*60);
                        Spike_TS_adj{i_trial,i_channel,i_unit} = temp_TS(trialidx{i_trial}) - 60*(i_trial-1);
                    end
                    
                end
            end

            % Save the Data to a .mat file
            fprintf('\n')
            fprintf('Saving: %s',  save_file);
            fprintf('\n')
            save(save_file,'LFP_TS_adj','LFP_mV_adj','Spike_TS_adj','Laser_TS_adj');
        end
    end
end