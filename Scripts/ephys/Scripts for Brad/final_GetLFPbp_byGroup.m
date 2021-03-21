% Get the LFP band power for each mouse at stimulation frequency and save them

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
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    freq = mouseList{i,4};
    
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
        
        % Calculate the band power
        fprintf('Calculating Band Power...\n');
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
                    
                    ts = LFP_mV_adj{i_trial,i_channel};
                    bp(i_trial,1) = bandpower(ts(b1),1000,[freq-0.01 freq+0.01]);
                    bp(i_trial,2) = bandpower(ts(b2),1000,[freq-0.01 freq+0.01]);
                    bp(i_trial,3) = bandpower(ts(b3),1000,[freq-0.01 freq+0.01]);
                    bp(i_trial,4) = bandpower(ts(b1),1000,[0 500]);
                    bp(i_trial,5) = bandpower(ts(b2),1000,[0 500]);
                    bp(i_trial,6) = bandpower(ts(b3),1000,[0 500]);
                end
                
                if i_channel <= 16
                    i_mouse = find(cellfun(@(x) isequal([mouse region],x),mouseTracker{i_group,1}),1);
                    if isempty(i_mouse)
                        mouseTracker{i_group,1}{end+1} = [mouse region];
                        i_mouse = size(mouseTracker{i_group,1},2);
                    end
                    All_Data{i_group,1,i_frequency,i_channel,i_mouse} = bp;
                else
                    i_mouse = find(cellfun(@(x) isequal(mouse,x),mouseTracker{i_group,i_region}),1);
                    if isempty(i_mouse)
                        mouseTracker{i_group,i_region}{end+1} = mouse;
                        i_mouse = size(mouseTracker{i_group,i_region},2);
                    end
                    All_Data{i_group,i_region,i_frequency,i_channel-16,i_mouse} = bp;
                end
                clear bp;
            end
        end
    end
end

fprintf('Saving: LFPbp.mat\n\n');
SaveFilename = fullfile(working_dir,'LFPbp.mat');
save(SaveFilename,'All_Data');