function final_RunQC2_BE(storage,base_fold,slash)
% Perform quality control step 2 on the rearranged data


%% 
stim = {'0_1' '0_5' '1_0'};

for i_mouse = 1:size(base_fold,1)
    QC_Status = [];
    
    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '_adj.mat'];
        save_file = [data_fold stim{i_stim} '_QC_2.txt'];
        clipidx_file = [data_fold stim{i_stim} '_clipidx.mat'];
        
        if exist(data_file,'file')
            load(data_file);
            msg = sprintf('Loading: %s\n',data_file);
            QC_Status = [QC_Status msg];
            
            % Check the laser timing
            msg = sprintf('   Checking laser timing\n');
            fprintf(msg);
            QC_Status = [QC_Status msg];
            laserFreq = 10;
            expectedLaser = 24:1/laserFreq:36;
            expectedLaser = expectedLaser(1:end-1);
            k_trial = 0;
            for i_trial = 1:length(Laser_TS_adj)
                if ~isequal(round(expectedLaser,3)',round(Laser_TS_adj{i_trial},3))
                    msg = sprintf('      Laser Error: Trial %d\n',i_trial);
                    fprintf(msg);
                    QC_Status = [QC_Status msg];
                else
                    k_trial = k_trial+1;
                end
            end
            
            if k_trial == length(Laser_TS_adj)
                msg = sprintf('      PASS\n');
                fprintf(msg);
                QC_Status = [QC_Status msg];
            end
                
            % Check for LFP timing fragments
            msg = sprintf('   Checking LFP timing for dropped fragments\n');
            fprintf(msg);
            QC_Status = [QC_Status msg];
            samplingFreq = 1000;
            expectedLength = samplingFreq * 60;
            k_trial = 0;
            for i_trial = 1:length(LFP_TS_adj)
                if ~isequal(expectedLength,length(LFP_TS_adj{i_trial}))
                    msg = sprintf('      Fragment Error: Trial %d\n',i_trial);
                    fprintf(msg);
                    QC_Status = [QC_Status msg];
                else
                    k_trial = k_trial+1;
                end
            end
            
            if k_trial == length(Laser_TS_adj)
                msg = sprintf('      PASS\n');
                fprintf(msg);
                QC_Status = [QC_Status msg];
            end
                
            % Check for clipping in the LFP
            msg = sprintf('   Checking for clipping\n');
            fprintf(msg);
            QC_Status = [QC_Status msg];
            k_mix = 0;
            clipidx = [];
            for i_trial = 1:size(LFP_mV_adj,1)
                for i_channel = 1:size(LFP_mV_adj,2)
                    temp_TS = LFP_mV_adj{i_trial,i_channel};
                    % Find consecutive max/min values
                    max_idx = find(temp_TS == max(temp_TS)); % Indices with max values
                    min_idx = find(temp_TS == min(temp_TS));
                    max_idx_diff = diff(max_idx); % Distance between indices
                    min_idx_diff = diff(min_idx);
                    max_clip = sum(max_idx_diff == 1); % Total number of adjacent indices with max values
                    min_clip = sum(min_idx_diff == 1);
                    if max_clip > 0 || min_clip > 0
                        msg = sprintf('      Clipping Error: Trial %d Channel %d >> %d max, %d min\n',i_trial,i_channel,max_clip,min_clip);
                        fprintf(msg);
                        QC_Status = [QC_Status msg];
                        clipidx = [clipidx; i_trial i_channel];
                    else
                        k_mix = k_mix + 1;
                    end
                end
            end
            
            if k_mix == size(LFP_mV_adj,1) * size(LFP_mV_adj,2)
                msg = sprintf('      PASS\n');
                fprintf(msg);
                QC_Status = [QC_Status msg];
            end
                
            % Save the results
            QC_Status = [QC_Status newline newline];
            fprintf('\n')
            fprintf('Saving: %s',  save_file);
            fprintf('\n')
            dlmwrite(save_file,QC_Status,'delimiter','');
            save(clipidx_file,'clipidx');
        end
    end
end