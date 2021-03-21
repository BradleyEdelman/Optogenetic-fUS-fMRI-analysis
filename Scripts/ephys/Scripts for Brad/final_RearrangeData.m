% Load and crop the data into 20 trials of 20" on, 20" on, 20" off
% Also compare the identified units across different recordings and
% frequencies, make sure you include them all

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

Mice = {'072419_A'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'6','10','20'};
Recordings = {'A','B','C','D'};
overwrite = 0; % 0 to just process new data, 1 for all data

%%
for i_mouse = 1:length(Mice)
    mouse = Mice{i_mouse};
    
    for i_region = 1:length(Regions)
        region = sprintf('BF_%s',Regions{i_region});
                
        for i_freq = 1:length(Frequencies)
            save_flag = 0;
            frequency = Frequencies{i_freq};
            
            for i_recording = 1:length(Recordings)
                recording = Recordings{i_recording};
                Filename = fullfile(working_dir,mouse,sprintf('%s_%sHz_%s.mat',region,frequency,recording));
                SaveFilename = fullfile(working_dir,mouse,sprintf('%s_%sHz_adj.mat',region,frequency));
                if exist(Filename,'file')
                    if overwrite == 1 || ~exist(SaveFilename,'file')
                        save_flag = 1;
                        fprintf('Loading: %s\n',Filename);
                        load(Filename);
                        
                        i_trial = 1+(i_recording-1)*5;
                        TS_delta = Laser_TS(1) - 20;
                        
                        % Create a cell array to hold the shifted Laser_TS values for each trial
                        temp_TS = Laser_TS - TS_delta;
                        trial_1_idx = find(temp_TS >= 0 & temp_TS < 60);
                        trial_2_idx = find(temp_TS >= 60 & temp_TS < 120);
                        trial_3_idx = find(temp_TS >= 120 & temp_TS < 180);
                        trial_4_idx = find(temp_TS >= 180 & temp_TS < 240);
                        trial_5_idx = find(temp_TS >= 240 & temp_TS < 300);
                        Laser_TS_adj{i_trial} = temp_TS(trial_1_idx);
                        Laser_TS_adj{i_trial+1} = temp_TS(trial_2_idx) - 60;
                        Laser_TS_adj{i_trial+2} = temp_TS(trial_3_idx) - 120;
                        Laser_TS_adj{i_trial+3} = temp_TS(trial_4_idx) - 180;
                        Laser_TS_adj{i_trial+4} = temp_TS(trial_5_idx) - 240;                       
                        
                        % Create a cell array to hold the shifted LFP_TS values for each trial
                        temp_TS = LFP_TS - TS_delta;
                        trial_1_idx = find(temp_TS >= 0 & temp_TS < 60);
                        trial_2_idx = find(temp_TS >= 60 & temp_TS < 120);
                        trial_3_idx = find(temp_TS >= 120 & temp_TS < 180);
                        trial_4_idx = find(temp_TS >= 180 & temp_TS < 240);
                        trial_5_idx = find(temp_TS >= 240 & temp_TS < 300);
                        LFP_TS_adj{i_trial} = temp_TS(trial_1_idx);
                        LFP_TS_adj{i_trial+1} = temp_TS(trial_2_idx) - 60;
                        LFP_TS_adj{i_trial+2} = temp_TS(trial_3_idx) - 120;
                        LFP_TS_adj{i_trial+3} = temp_TS(trial_4_idx) - 180;
                        LFP_TS_adj{i_trial+4} = temp_TS(trial_5_idx) - 240;
                        % Create a cell array for LFP_mV aligned to the laser with all trials
                        for i_channel = 1:length(LFP_mV)
                            temp_mV = LFP_mV{i_channel};
                            LFP_mV_adj{i_trial,i_channel} = temp_mV(trial_1_idx);
                            LFP_mV_adj{i_trial+1,i_channel} = temp_mV(trial_2_idx);
                            LFP_mV_adj{i_trial+2,i_channel} = temp_mV(trial_3_idx);
                            LFP_mV_adj{i_trial+3,i_channel} = temp_mV(trial_4_idx);
                            LFP_mV_adj{i_trial+4,i_channel} = temp_mV(trial_5_idx);
                        end
                        
                        % Create a Trial x Channel x Unit cell array for Spike_TS aligned to the laser
                        for i_channel = 1:size(Spike_TS,2)
                            for i_unit = 1:size(Spike_TS,3)
                                temp_TS = Spike_TS{1,i_channel,i_unit} - TS_delta;
                                trial_1_idx = find(temp_TS >= 0 & temp_TS < 60);
                                trial_2_idx = find(temp_TS >= 60 & temp_TS < 120);
                                trial_3_idx = find(temp_TS >= 120 & temp_TS < 180);
                                trial_4_idx = find(temp_TS >= 180 & temp_TS < 240);
                                trial_5_idx = find(temp_TS >= 240 & temp_TS < 300);
                                Spike_TS_adj{i_trial,i_channel,i_unit} = temp_TS(trial_1_idx);
                                Spike_TS_adj{i_trial+1,i_channel,i_unit} = temp_TS(trial_2_idx) - 60;
                                Spike_TS_adj{i_trial+2,i_channel,i_unit} = temp_TS(trial_3_idx) - 120;
                                Spike_TS_adj{i_trial+3,i_channel,i_unit} = temp_TS(trial_4_idx) - 180;
                                Spike_TS_adj{i_trial+4,i_channel,i_unit} = temp_TS(trial_5_idx) - 240;
                            end
                        end
                    end
                end
            end

            % Save the Data to a .mat file
            if save_flag == 1
                fprintf('Saving: %s  %s_%sHz_adj.mat\n\n',mouse,region,frequency);
                save(SaveFilename,'LFP_TS_adj','LFP_mV_adj','Spike_TS_adj','Laser_TS_adj');
            end
        end
    end
end