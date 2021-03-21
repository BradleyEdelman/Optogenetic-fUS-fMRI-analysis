function final_CheckSpikesDuringClips_BE(storage,base_fold,slash)
% Investigate how clipping affects spiking behavior
fs = 1000; % sampling frequency


%%
stim = {'0_1' '0_5' '1_0'};

for i_mouse = 1:size(base_fold,1)
    
    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '.mat'];
%         save_file = [data_fold stim{i_stim} '_QC_2.mat'];
        
        if exist(data_file,'file')
            load(data_file);
                    
            for i_ch = 1:32
                % Find possible clipping
                mV = LFP_mV{i_ch};
                max_mV = max(mV);
                min_mV = min(mV);
                lim_mV = max(max_mV,abs(min_mV));
                max_mV = lim_mV;
                min_mV = -lim_mV;
                max_idx = find(mV == max_mV);
                min_idx = find(mV == min_mV);
                if (length(max_idx)+length(min_idx)) == 1
                    max_idx = [];
                    min_idx = [];
                end
                max_pts = LFP_TS(max_idx);
                min_pts = LFP_TS(min_idx);

                clear pts diff_pts diff_idx;
                pts{1} = max_pts;
                pts{2} = min_pts;

                diff_pts{1} = diff(floor(max_pts*1000));
                diff_pts{2} = diff(floor(min_pts*1000));
                diff_idx{1} = find(diff_pts{1} ~= 1) + 1;
                diff_idx{2} = find(diff_pts{2} ~= 1) + 1;
                d_list{1} = [];
                d_list{2} = [];
                final_list = [];
                
                for i=1:2
                    d_start = [];
                    d_end = [];
                    if ~isempty(pts{i})
                        d_start = pts{i}(1);
                        if isempty(diff_pts{i})
                            d_end = d_start;
                        elseif isempty(diff_idx{i})
                            d_end = pts{i}(end);
                        else
                            d_end = pts{i}(diff_idx{i}(1)-1);
                            if (length(diff_idx{i}) > 1)
                                for i_diff = 1:(length(diff_idx{i})-1)
                                    d_start = [d_start  pts{i}(diff_idx{i}(i_diff))];
                                    d_end = [d_end  pts{i}(diff_idx{i}(i_diff+1)-1)];
                                end
                            end
                            d_start = [d_start  pts{i}(diff_idx{i}(end))];
                            d_end = [d_end  pts{i}(end)];
                        end
                    end
                    for k = 1:length(d_start)
                        d_list{i} = [d_list{i};[d_start(k),d_end(k)]];
                    end
                end
                
                final_list = [d_list{1};d_list{2}];
                if ~isempty(final_list)
                    final_list(:,1) = final_list(:,1) - 1/fs; % to account for spikes between LFP_TS points
                    final_list(:,2) = final_list(:,2) + 1/fs; % to account for spikes between LFP_TS points
                end
                
                for i_unit = 1:6
                    temp_TS = Spike_TS{1,i_ch,i_unit};
                    num_spikes = length(temp_TS);
                    counter = 0;
                    if any(ismember([pts{1} pts{2}],Spike_TS{1,i_ch,i_unit}))
                        fprintf('   :[%d,%d]\n',i_ch,i_unit);
                    end
                    
                    for i_fl = 1:size(final_list,1)
                        for i_TS = 1:length(temp_TS)
                            if ((temp_TS(i_TS) >= final_list(i_fl,1)) && (temp_TS(i_TS) <= final_list(i_fl,2)))
                                counter = counter + 1;
                                if ((counter/num_spikes*100) >= 1)
                                    fprintf('   %s:[%d,%d]#%d - %.2f%%\n',recording,i_ch,i_unit-1,i_TS,counter/num_spikes*100);
                                end
                            end
                        end
                    end
                end
                
            end
                    
        end
    end
end
close all;