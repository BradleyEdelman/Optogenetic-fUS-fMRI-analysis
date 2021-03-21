function final_PlotLFPspectrograms_Averaged_BE(storage,base_fold,slash)
% Plot and save the LFP Spectrogram averaged across trials


t_window = 1000;
t_overlap = 800;
nfft = 500;
fs = 1000;
climit = [-75 -35];
ylimit = [0 80];

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
Fig2 = figure('Position', get(0, 'Screensize'),'visible','off');
stim = {'0_1' '0_5' '1_0'};

for i_mouse = 1:size(base_fold,1)

    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '_clean.mat'];
    
        if exist(data_file,'file')
            load(data_file);
        
            n_channels = size(LFP_mV_adj,2);
            n_trials = size(LFP_mV_adj,1);
        
            % Calculate the averages
            fprintf('Averaging LFP...\n');
            for i_channel = 1:n_channels
                if ~isempty(vertcat(LFP_mV_adj{:,i_channel}))
                    for i_trial = 1:n_trials
                        if isempty(LFP_mV_adj{i_trial,i_channel})
                        elseif length(LFP_mV_adj{i_trial,i_channel}) > 60000
                            LFP_mV_adj{i_trial,i_channel} = LFP_mV_adj{i_trial,i_channel}(1:60000);
                        elseif length(LFP_mV_adj{i_trial,i_channel}) < 60000
                            endpoint = length(LFP_mV_adj{i_trial,i_channel});
                            LFP_mV_adj{i_trial,i_channel}((endpoint+1):60000) = LFP_mV_adj{i_trial,i_channel}(endpoint);
                        end
                    end
                    trials_per_channel = horzcat(LFP_mV_adj{:,i_channel})';
                    
                    if n_trials > 10
                        trials_per_channel(1,:) = [];
                    end
                    
                    average_trials{1,i_channel} = mean(trials_per_channel,1);
                else
                    average_trials{1,i_channel} = nan(1,60000);
                end
            end
        
            % Plot the spectrograms of the averages per channel
            fprintf('Plotting/Saving...\n');
            for i_channel = 1:n_channels
                figure(Fig1); subplot(4,8,i_channel);
                if sum(isnan(average_trials{1,i_channel})) == 0
                    spectrogram(average_trials{1,i_channel},t_window,t_overlap,nfft,fs,'yaxis');
                    hold on;
                    text(0.02,0.90,['Ch.' char(num2str(i_channel))],'Units','Normalized','FontWeight','bold','Color','w');
                    % Format the plots
                    ax=gca; ax.CLim = climit;
                    ylim(ylimit);
                end
            end
            % Save the plot
            saveas(Fig1,[data_fold 'PlotLFPspectrograms_Averaged_' stim{i_stim} '.png']);
            clf;
            
            % Plot the spectrograms of the averages per region
            average_trials_region{1} = nanmean(vertcat(average_trials{1,1:8}),1);
            average_trials_region{2} = nanmean(vertcat(average_trials{1,9:16}),1);
            average_trials_region{3} = nanmean(vertcat(average_trials{1,17:24}),1);
            average_trials_region{4} = nanmean(vertcat(average_trials{1,25:32}),1);
            t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
            for i_region = 1:4
                figure(Fig2); subplot(2,2,i_region);
                spectrogram(average_trials_region{i_region},t_window,t_overlap,nfft,fs,'yaxis');
                hold on;
                text(0.02,0.90,[t_bank{i_region}],'Units','Normalized','FontWeight','bold','Color','w');
                % Format the plots
                ax=gca; ax.CLim = climit;
                ylim(ylimit);
            end
            % Save the plot
            saveas(Fig1,[data_fold 'PlotLFPspectrograms_Averaged_Region_' stim{i_stim} '.png']);
            clf;
            
        end
    end    
end