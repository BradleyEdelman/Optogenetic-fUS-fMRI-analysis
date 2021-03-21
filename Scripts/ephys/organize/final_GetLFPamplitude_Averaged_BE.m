function final_GetLFPamplitude_Averaged_BE(storage,base_fold,slash)
% Get the LFP amplitude for each mouse at stimulation frequency and save them


%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
Fig2 = figure('Position', get(0, 'Screensize'),'visible','off');

% 24 second baseline, 12 second stim, 24 second baseline
b1 = 1:24000;
b2 = 24001:36000;
b3 = 36001:60000;

fs = 1000;
ftype = 'bandpass';

stim = {'0_1' '0_5' '1_0'};

for i_mouse = 1:size(base_fold,1)
    
    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '_clean.mat'];
        save_file = [data_fold stim{i_stim} '_LFP_amp.mat'];
        
        if exist(data_file,'file')
            load(data_file);
    
            fc1 = 8; % highpass frequency
            fc2 = 12; % lowpass frequency
            Wn = [fc1 fc2]*2/fs;
            [b,a] = butter(4, Wn, ftype);
        
            n_channels = size(LFP_mV_adj,2);
            n_trials = size(LFP_mV_adj,1);
                
            % Calculate the LFP amplitude
            fprintf('  Calculating LFP amplitude...\n');
            for i_channel = 1:n_channels
                if ~isempty(LFP_mV_adj(:,i_channel))

                    for i_trial = 1:n_trials
                        if isempty(LFP_mV_adj{i_trial,i_channel})
                        elseif length(LFP_mV_adj{i_trial,i_channel}) > 60000
                            LFP_mV_adj{i_trial,i_channel} = LFP_mV_adj{i_trial,i_channel}(1:60000);
                        elseif length(LFP_mV_adj{i_trial,i_channel}) < 60000
                            endpoint = length(LFP_mV_adj{i_trial,i_channel});
                            LFP_mV_adj{i_trial,i_channel}((endpoint+1):60000) = LFP_mV_adj{i_trial,i_channel}(endpoint);
                        end
                        
                        if ~isempty(LFP_mV_adj{i_trial,i_channel})
                            % Band pass filter the data
                            ts = LFP_mV_adj{i_trial,i_channel};
                            ts_filt = filtfilt(b, a, ts);

                            if any(isnan(ts_filt))
                                warning('Filter failure on Ch.%d Trial %d',i_channel,i_trial);
                            end

                            ts = abs(ts);
                            ts_filt = abs(ts_filt);

                            amp(i_trial,1) = mean(ts_filt(b1));
                            amp(i_trial,2) = mean(ts_filt(b2));
                            amp(i_trial,3) = mean(ts_filt(b3));
                            
                            mean_base{i_trial} = mean(ts_filt(b1));
                            std_base{i_trial} = std(ts_filt(b1));
                            Stim{i_trial} = ts_filt(b2);
                            Ave_z_trial{i_trial} = mean((Stim{i_trial} - mean_base{i_trial})/std_base{i_trial});
                            
                            amp(i_trial,4) = mean(ts(b1));
                            amp(i_trial,5) = mean(ts(b2));
                            amp(i_trial,6) = mean(ts(b3));
                        else
                            amp(i_trial,1:6) = nan;
                        end
                    end

                    LFP_amp{i_channel} = amp;
                    
                    % percent baseline change
                    amp = LFP_amp{i_channel};
                    amp_change = amp(:,2)./amp(:,1)*100 - 100;
                    amp_change = (amp(:,2) - amp(:,1))./amp(:,1)*100;
                    avg = mean(amp_change);
                    sd = std(amp_change);
                    sem = sd / sqrt(size(amp_change,1));
                    LFP_amp_change(1,i_channel) = avg;
                    LFP_amp_change(2,i_channel) = sem;
                    
                    LFP_ave_z(1,i_channel) = mean(cell2mat(Ave_z_trial));
                    clear amp
                    
                end
            end
            
            % For each of the four regions recorded from
            %('L CPu' 'L Cortex' 'R CPu' 'R Cortex')
            r_idx = {1:8, 9:16, 17:24, 25:32};
            for i_region = 1:4
                for i_trial = 1:size(LFP_mV_adj,1)
                    if ~isempty(horzcat(LFP_mV_adj{i_trial,r_idx{i_region}}))
                        ts_region = horzcat(LFP_mV_adj{i_trial,r_idx{i_region}});
                        ts_region = mean(ts_region,2);

                        ts_filt = filtfilt(b, a, ts_region);
                        if any(isnan(ts_filt))
                            warning('Filter failure on Region %d',i_region);
                        end
                        ts = abs(ts_region);
                        ts_filt = abs(ts_filt);

                        amp_r(i_trial,1) = nanmean(ts_filt(b1));
                        amp_r(i_trial,2) = nanmean(ts_filt(b2));
                        amp_r(i_trial,3) = nanmean(ts_filt(b3));
                        
%                         if i_stim == 2
%                             x = 2;
%                         end
                        
                        mean_base_r{i_trial} = mean(ts_filt(b1));
                        std_base_r{i_trial} = std(ts_filt(b1));
                        Stim_r{i_trial} = ts_filt(b2);
                        Ave_z_trial_r{i_trial} = mean((Stim_r{i_trial} - mean_base_r{i_trial})/std_base_r{i_trial});
                        
                        amp_r(i_trial,4) = nanmean(ts(b1));
                        amp_r(i_trial,5) = nanmean(ts(b2));
                        amp_r(i_trial,6) = nanmean(ts(b3));
                    else
                        amp_r(i_trial,1:6) = nan;
                    end
                end
                
                LFP_amp_region{i_region} = amp_r;
                
                % percent baseline change
                amp = LFP_amp_region{i_region};
                amp_change = amp(:,2)./amp(:,1)*100 - 100;
                amp_change = (amp(:,2) - amp(:,1))./amp(:,1)*100;
                avg = nanmean(amp_change);
                sd = nanstd(amp_change);
                sem = sd / sqrt(size(amp_change,1));
                LFP_amp_region_change(1,i_region) = avg;
                LFP_amp_region_change(2,i_region) = sem;
                
                LFP_ave_z_region(1,i_region) = mean(cell2mat(Ave_z_trial_r));
                clear amp
            end
            
            fprintf('Saving: %s\n\n',save_file);
            save(save_file,'LFP_amp','LFP_amp_change','LFP_ave_z','LFP_amp_region','LFP_amp_region_change','LFP_ave_z_region');
            clear amp
            
        end
    end
end

