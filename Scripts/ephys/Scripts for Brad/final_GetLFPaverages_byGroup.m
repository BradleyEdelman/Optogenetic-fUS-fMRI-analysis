% Get the LFP averages for each mouse and save them

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

fileTag = '_filt';

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');

All_Data = cell(length(Groups),length(Regions),2,16);
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
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
                average_trials{1,i_channel} = mean(trials_per_channel,1);
                
                if i_channel <= 16
                    All_Data{i_group,1,i_frequency,i_channel} = [All_Data{i_group,1,i_frequency,i_channel};average_trials{1,i_channel}];
                else
                    All_Data{i_group,i_region,i_frequency,i_channel-16} = [All_Data{i_group,i_region,i_frequency,i_channel-16};average_trials{1,i_channel}];
                end
                
            end
        end
    end
end

fprintf('Saving: LFPaverages.mat\n\n');
SaveFilename = fullfile(working_dir,'LFPaverages.mat');
save(SaveFilename,'All_Data');