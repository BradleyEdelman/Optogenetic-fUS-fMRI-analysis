% Get the LFP PSD for each mouse and save them

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

fileTag = '_filt';

%%
t_window = 1000;
t_overlap = 800;
nfft = 1000;
fs = 1000;

All_Data = cell(length(Groups),length(Regions),2,16);
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
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
        
        % Calculate the averages
        fprintf('Computing/Averaging PSD...\n');
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
                    [s,f,t,ps] = spectrogram(LFP_mV_adj{i_trial,i_channel},t_window,t_overlap,nfft,fs,'yaxis','power');
                    total_ps = sum(sum(ps,2),1);
                    p(:,:,i_trial) = ps/total_ps;
                end
                average_trials{1,i_channel} = mean(p,3);
                clear p;
                
                if i_channel <= 16
                    All_Data{i_group,1,i_frequency,i_channel} = cat(3,All_Data{i_group,1,i_frequency,i_channel},average_trials{1,i_channel});
                else
                    All_Data{i_group,i_region,i_frequency,i_channel-16} = cat(3,All_Data{i_group,i_region,i_frequency,i_channel-16},average_trials{1,i_channel});
                end
                
            end
        end
    end
end

for i_group = 1:length(Groups)
    fprintf('Saving: LFPpsd_%s.mat\n\n',Groups{i_group});
    SaveFilename = fullfile(working_dir,sprintf('LFPpsd_%s.mat',Groups{i_group}));
    Data = All_Data(i_group,:,:,:);
    save(SaveFilename,'Data');
end