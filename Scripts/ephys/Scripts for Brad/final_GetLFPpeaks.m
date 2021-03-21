% Get the delays between the laser onset and first LFP peak and save them

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

fileTag = '_filt';
range = 19501:20500;

%%
All_Data = cell(length(Groups),length(Regions),2,32);
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});

    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
    if exist(Filename,'file')
        fprintf('%d - Loading: %s\n',i,Filename);
        load(Filename);

        [~,i_group] = ismember(group,Groups);
        [~,i_region] = ismember(region,Regions);
        if frequency == '20'
            i_frequency = 2;
        else
            i_frequency = 1;
        end

        num_trials = size(LFP_mV_adj,1);
        num_channels = size(LFP_mV_adj,2);

        % Loop through
        for i_ch = 1:num_channels
            
            for xx = 1:size(LFP_mV_adj,1)
                if isempty(LFP_mV_adj{xx,i_ch})
                    LFP_mV_adj{xx,i_ch} = zeros(1,60000);
                end
            end

            % Filter the data
            data = cellfun(@(x) x(range),LFP_mV_adj(:,i_ch),'UniformOutput',false);
            data = cell2mat(data);
            data = reshape(data,[],num_trials)';
            data = mean(data,1);

            All_Data{i_group,i_region,i_frequency,i_ch} = [All_Data{i_group,i_region,i_frequency,i_ch};data];

        end

    end
end

fprintf('Saving: LFPpeaks.mat\n\n');
SaveFilename = fullfile(working_dir,'LFPpeaks.mat');
save(SaveFilename,'All_Data','range');