% Perform quality control step 2 on the rearranged data

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

Mice = {'072419_A'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'6','10','20'};

%%
for i_mouse = 1:length(Mice)
    mouse = Mice{i_mouse};
    for i_region = 1:length(Regions)
        region = sprintf('BF_%s',Regions{i_region});
        for i_freq = 1:length(Frequencies)
            frequency = Frequencies{i_freq};
            QC_Status = [];
            Filename = fullfile(working_dir,mouse,sprintf('%s_%sHz_adj.mat',region,frequency));
            if exist(Filename,'file')
                msg = sprintf('Loading: %s\n',Filename);
                fprintf('Loading: %s\n',Filename);
                load(Filename);
                QC_Status = [QC_Status msg];
                
                % Check the laser timing
                msg = sprintf('   Checking laser timing\n');
                fprintf(msg);
                QC_Status = [QC_Status msg];
                laserFreq = str2double(frequency);
                expectedLaser = 20:1/laserFreq:40;
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
                msg = fprintf('\n\n');
                QC_Status = [QC_Status newline newline];
                OutputFilename = fullfile(working_dir,mouse,sprintf('QC2_%s_%sHz.txt',region,frequency));
                dlmwrite(OutputFilename,QC_Status,'delimiter','');
            end
        end
    end
end