% Investigate how clipping affects spiking behavior

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';
Mice = {'061519_A','061519_B','061519_C','061619_C','061719_A','061719_B','061719_C','061819_A','061819_B','061819_C','061819_D','062019_A','062019_B','062019_C','062019_D','062119_A','062119_B','063019_A','063019_B','063019_C','070119_A','070119_B','070119_C','070619_A','070619_B','070619_C','070719_A','070719_B','070719_C','070719_D','070819_A','070819_B','070819_C'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'6','10','20'};
Recordings = {'A','B','C','D'};
fs = 1000; % sampling frequency


Mice = {'072419_A'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'10','20'};
Recordings = {'A','B','C','D'};


%%
for i_mouse = 1:length(Mice)
    mouse = Mice{i_mouse};
    for i_region = 1:length(Regions)
        region = sprintf('BF_%s',Regions{i_region});
        for i_freq = 1:length(Frequencies)
            frequency = Frequencies{i_freq};
            for i_recording = 1:length(Recordings)
                recording = Recordings{i_recording};
                Filename = fullfile(working_dir,mouse,sprintf('%s_%sHz_%s.mat',region,frequency,recording));
                if exist(Filename,'file')
                    fprintf('Loading: %s\n',Filename);
                    load(Filename);
                    
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
                                fprintf('   %s:[%d,%d]\n',recording,i_ch,i_unit);
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
    end
end
close all;