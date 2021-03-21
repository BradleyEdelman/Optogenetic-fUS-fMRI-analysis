% Get the LFP averages for each mouse and save them


%%
Fig1 = figure;%('Position', get(0, 'Screensize'),'visible','off');

stim = {'0_1' '0_5' '1_0'};

for i_stim = 1:size(stim,2)
    save_fold = [storage stim{i_stim} slash];
    if ~exist(save_fold,'dir'); mkdir(save_fold); end
    save_file_channel = [save_fold 'Group_channel.mat'];
    save_file_region = [save_fold 'Group_region.mat'];
    
    for i_mouse = 1:size(base_fold,1)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '.mat'];
        
        if exist(data_file,'file')
            load(data_file)
        
            n_channels = size(LFP_mV_adj,2);
            n_trials = size(LFP_mV_adj,1);
            
            % Calculate the averages
            fprintf('Averaging LFP...\n');
            for i_channel = 1:n_channels
                if ~isempty(LFP_mV_adj{1,i_channel})
                    for i_trial = 1:n_trials
                        if length(LFP_mV_adj{i_trial,i_channel}) > 60000
                            LFP_mV_adj{i_trial,i_channel} = LFP_mV_adj{i_trial,i_channel}(1:60000);
                        end
                        if length(LFP_mV_adj{i_trial,i_channel}) < 60000
                            endpoint = length(LFP_mV_adj{i_trial,i_channel});
                            LFP_mV_adj{i_trial,i_channel}((endpoint+1):60000) = LFP_mV_adj{i_trial,i_channel}(endpoint);
                        end
                    end
                    trials_per_channel = reshape(cell2mat(LFP_mV_adj(:,i_channel)),length(LFP_mV_adj{1,1}),n_trials)';
                    average_trials{i_mouse,i_channel} = mean(trials_per_channel,1);
                end
            end 
        end
    end
    
    average_trials_channel = mean(vertcat(average_trials
    
end

fprintf('Saving: LFPaverages.mat\n\n');
SaveFilename = fullfile(working_dir,'LFPaverages.mat');
save(SaveFilename,'All_Data');