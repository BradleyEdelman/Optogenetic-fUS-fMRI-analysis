% Count the delays between the laser onset and the spikes

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

fileTag = '_filt';

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
        
        % Calculate the spike delays
        Laser_delta = Laser_TS_adj{1}(2) - Laser_TS_adj{1}(1);
        Spike_Delays = cellfun(@(x) mod(x(x>=20 & x<=40),Laser_delta),Spike_TS_adj,'UniformOutput',false);
        with_spikes_idx = ~cellfun('isempty',Spike_Delays);
        First_Delays = cell(size(Spike_Delays));
        First_Delays(with_spikes_idx) = cellfun(@(x) x(1),Spike_Delays(with_spikes_idx),'UniformOutput',false);
        First_Delays(~with_spikes_idx) = cellfun(@(x) 0,First_Delays(~with_spikes_idx),'UniformOutput',false);
        First_Delays = cell2mat(First_Delays);
        
        % Save the results
        fprintf('Saving: %s  BF_%s_%sHz_SpikeDelays.mat\n\n',mouse,region,frequency);
        SaveFilename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_SpikeDelays.mat',region,frequency));
        save(SaveFilename,'Spike_Delays','First_Delays');
    end
end