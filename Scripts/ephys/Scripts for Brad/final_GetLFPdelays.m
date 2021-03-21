% Count the delays between the laser onset and first LFP peak

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

fileTag = '_filt';
maxThresh = [60,80,90,95,100]; % Percent of max height to look for

%%
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        
        % Calculate the LFP delays
        fprintf('Calculating LFP Delays...\n');
        LFP_Delays = cell(size(LFP_mV_adj,1),size(LFP_mV_adj,2),length(maxThresh));
        for i_channel = 1:size(LFP_mV_adj,2)
            if ~isempty(LFP_mV_adj{1,i_channel})
                for i_trial = 1:size(LFP_mV_adj,1)
                    mV = LFP_mV_adj{i_trial,i_channel};
                    mV = abs(mV); % Just look at absolute value
                    endpoint = 20000+floor(1000/str2num(frequency));
                    mV = mV(20000:endpoint); % Just look at the first stimulation
                    for i_thresh = 1:length(maxThresh)
                        thresh = maxThresh(i_thresh);
                        thresh = max(mV)*thresh/100;
                        idx = find(mV >= thresh,1);
                        LFP_Delays{i_trial,i_channel,i_thresh} = idx - 1;
                    end
                end
            end
        end
        
        LFP_Delays(cellfun('isempty',LFP_Delays)) = cellfun(@(x) 0,LFP_Delays(cellfun('isempty',LFP_Delays)),'UniformOutput',false);
        LFP_Delays = cell2mat(LFP_Delays);
        
        % Save the results
        fprintf('Saving: %s  BF_%s_%sHz_LFPdelays.mat\n\n',mouse,region,frequency);
        SaveFilename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_LFPdelays.mat',region,frequency));
        save(SaveFilename,'LFP_Delays','maxThresh');
    end
end