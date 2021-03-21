% Define bad channels/trials/units to filter out of the data

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';
[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');

rows = 2:137;
drop_unsorted = 1; % To drop the 'unsorted' units, located at (:,:,1)

%%

% Load data from xlsx file
mouseList = raw(rows,1:4);
mouseList(:,5) = cellfun(@(x) double(strsplit(string(x),',')),raw(rows,21),'UniformOutput',false);
mouseList(:,6) = cellfun(@(x) double(strsplit(string(x),',')),raw(rows,22),'UniformOutput',false);
mouseList(:,7) = cellfun(@(x) strsplit(string(x),';'),raw(rows,23),'UniformOutput',false);

mouseList(:,8) = raw(rows,5);
clear raw;

for i = 1:size(mouseList,1)
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    freq = num2str(mouseList{i,4}); 
    adj_Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_adj.mat',region,freq));
    filt_Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_filt.mat',region,freq));
    
    % Load the data
    fprintf('Loading: %s\n',adj_Filename);
    load(adj_Filename);
    
    if isequal(mouseList{i,8},'M180')
        % Reverse the channel order for the Poly2 electrode (M180)
        fprintf('   Reversing Poly2 Electrode Channels\n');
        LFP_temp = cell(size(LFP_mV_adj));
        Spike_temp = cell(size(Spike_TS_adj));
        for i_ch = 1:16
            LFP_temp(:,i_ch) = LFP_mV_adj(:,17-i_ch);
            Spike_temp(:,i_ch,:) = Spike_TS_adj(:,17-i_ch,:);
        end
        LFP_temp(:,17:32) = LFP_mV_adj(:,17:32);
        Spike_temp(:,17:32,:) = Spike_TS_adj(:,17:32,:);
        LFP_mV_adj = LFP_temp;
        Spike_TS_adj = Spike_temp;
    elseif isequal(mouseList{i,8},'M452')
        % Re-space the sites on M452 to match M180
        fprintf('   Re-spacing M452 Electrode Channels\n');
        LFP_temp = cell(size(LFP_mV_adj));
        Spike_temp = cell(size(Spike_TS_adj));
        tt = 1;
        for M180_ch = [1 3 5 7 9 11 13 15 16]
            M452_ch = 1:9;
            LFP_temp(:,M180_ch) = LFP_mV_adj(:,M452_ch(tt));
            Spike_temp(:,M180_ch,:) = Spike_TS_adj(:,M452_ch(tt),:);
            tt = tt+1;
        end
        LFP_temp(:,17:32) = LFP_mV_adj(:,17:32);
        Spike_temp(:,17:32,:) = Spike_TS_adj(:,17:32,:);
        LFP_mV_adj = LFP_temp;
        Spike_TS_adj = Spike_temp;
    elseif isequal(mouseList{i,8},'L117')
        % Re-space the sites on L117 to match M180
        fprintf('   Re-spacing L117 Electrode Channels\n');
        LFP_temp = cell(size(LFP_mV_adj));
        Spike_temp = cell(size(Spike_TS_adj));
        tt = 1;
        for M180_ch = [1 7 13 16]
            L117_ch = 1:4;
            LFP_temp(:,M180_ch) = LFP_mV_adj(:,L117_ch(tt));
            Spike_temp(:,M180_ch,:) = Spike_TS_adj(:,L117_ch(tt),:);
            tt = tt+1;
        end
        LFP_temp(:,17:32) = LFP_mV_adj(:,17:32);
        Spike_temp(:,17:32,:) = Spike_TS_adj(:,17:32,:);
        LFP_mV_adj = LFP_temp;
        Spike_TS_adj = Spike_temp;
    end

    % Filter the data
    clear badChannels badTrials badUnits
    if ~isnan(mouseList{i,5})
        fprintf('   Dropping Channels: ');
        badChannels = mouseList{i,5};
        for k = length(badChannels):-1:1
            bC = badChannels(k);
            fprintf('%d ',bC);
            for i_trial = 1:20
                LFP_mV_adj{i_trial,bC} = [];
                for i_unit = 1:size(Spike_TS_adj,3)
                    Spike_TS_adj{i_trial,bC,i_unit} = [];
                end
            end
        end
        fprintf(newline);
    else
        badChannels = [];
    end
    
    if ~isnan(mouseList{i,6})
        fprintf('   Dropping Trials: ');
        badTrials = mouseList{i,6};
        for k = length(badTrials):-1:1
            bT = badTrials(k);
            fprintf('%d ',bT);
            Laser_TS_adj(bT) = [];
            LFP_TS_adj(bT) = [];
            LFP_mV_adj(bT,:) = [];
            Spike_TS_adj(bT,:,:) = [];        
        end
        fprintf(newline);
    else
        badTrials = [];
    end
    
    d = 1;
    if drop_unsorted == 1
        fprintf('   Dropping Units: All Unsorted (unit 0)\n');
        Spike_TS_adj(:,:,1) = [];
        d = 0; % Compensate for all the unit numbers shifting down by one after deletion
    end
    if ~all(ismissing(mouseList{i,7}))
        fprintf('   Dropping Units: ');
        badUnits = mouseList{i,7};
        for k = length(badUnits):-1:1
            curr_unit = strsplit(badUnits(k),',');
            ch = char(curr_unit(1));
            ch = str2double(ch(2:end));
            unit = char(curr_unit(2));
            unit = str2double(unit(1:end-1));
            unit = unit + d;
            fprintf('[%d,%d] ',ch,unit);
            for zz = 1:size(Spike_TS_adj(:,ch,unit),1)
                Spike_TS_adj{zz,ch,unit} = [];
            end
        end
        fprintf(newline);
    else
        badUnits = [];
    end
    
    % Save the data to a new file with the filters saved as well
    fprintf('Saving: %s\n\n',filt_Filename);
    save(filt_Filename,'LFP_TS_adj','LFP_mV_adj','Spike_TS_adj','Laser_TS_adj','badChannels','badTrials','badUnits');
end