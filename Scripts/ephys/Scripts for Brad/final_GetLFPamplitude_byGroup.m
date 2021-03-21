% Get the LFP amplitude for each mouse at stimulation frequency and save them

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

fileTag = '_filt';

%%
All_Data = cell(length(Groups),length(Regions),2,16,2);
mouseTracker = cell(length(Groups),length(Regions));
for aa = 1:length(Groups)
    for bb = 1:length(Regions)
        mouseTracker{aa,bb} = {};
    end
end
b1 = 1:20000;
b2 = 20001:40000;
b3 = 40001:60000;

fs = 1000;
ftype = 'bandpass';

for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    freq = mouseList{i,4};
    
    fc1 = freq-2; % highpass frequency
    fc2 = freq+2; % lowpass frequency
    Wn = [fc1 fc2]*2/fs;
    [b,a] = butter(4, Wn, ftype);
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
    if exist(Filename,'file')
        fprintf('%d - Loading: %s\n',i,Filename);
        load(Filename);
        
        n_channels = size(LFP_mV_adj,2);
        n_trials = size(LFP_mV_adj,1);
        
        [~,i_group] = ismember(group,Groups);
        [~,i_region] = ismember(region,Regions);
        if frequency == '20'
            i_frequency = 2;
        else
            i_frequency = 1;
        end
                
        % Calculate the LFP amplitude
        fprintf('  Calculating LFP amplitude...\n');
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
                    amp(i_trial,4) = mean(ts(b1));
                    amp(i_trial,5) = mean(ts(b2));
                    amp(i_trial,6) = mean(ts(b3));
                end
                                
                if i_channel <= 16
                    i_mouse = find(cellfun(@(x) isequal([mouse region],x),mouseTracker{i_group,1}),1);
                    if isempty(i_mouse)
                        mouseTracker{i_group,1}{end+1} = [mouse region];
                        i_mouse = size(mouseTracker{i_group,1},2);
                    end
                    All_Data{i_group,1,i_frequency,i_channel,i_mouse} = amp;
                else
                    i_mouse = find(cellfun(@(x) isequal(mouse,x),mouseTracker{i_group,i_region}),1);
                    if isempty(i_mouse)
                        mouseTracker{i_group,i_region}{end+1} = mouse;
                        i_mouse = size(mouseTracker{i_group,i_region},2);
                    end
                    All_Data{i_group,i_region,i_frequency,i_channel-16,i_mouse} = amp;
                end
                clear amp;
            end
        end
    end
end

fprintf('Saving: LFPamplitude.mat\n\n');
SaveFilename = fullfile(working_dir,'LFPamplitude.mat');
save(SaveFilename,'All_Data');