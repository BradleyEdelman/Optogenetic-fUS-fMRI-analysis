% Count the SLOs in the LFP data

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

fileTag = '_filt';

%%

% % minT = 500; % Rest period before SLO (ms)
% % postT = 500; % Duriation of SLO measurement (ms)
% % postThresh = 0.02; % Threshold for SLO during postT
% % 
% % fc1 = 2; % highpass frequency
% % fc2 = 200; % lowpass frequency
% % fs = 1000;
% % Wn = [fc1 fc2]*2/fs;
% % ftype = 'bandpass';
% % [b,a] = butter(5, Wn, ftype);
% % 
% % for i = 1:size(mouseList,1)
% %     group = mouseList{i,1};
% %     mouse = mouseList{i,2};
% %     region = mouseList{i,3};
% %     frequency = num2str(mouseList{i,4});
% %     
% %     Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
% %     if exist(Filename,'file')
% %         fprintf('Loading: %s\n',Filename);
% %         load(Filename);
% %         
% %         SLOcount = zeros(size(LFP_mV_adj));
% %         SLOidx = cell(size(LFP_mV_adj));        
% %                 
% %         for xx = 1:size(LFP_mV_adj,1)
% %             for yy = 1:size(LFP_mV_adj,2)
% %                 if isempty(LFP_mV_adj{xx,yy})
% %                     LFP_mV_adj{xx,yy} = zeros(1,60000);
% %                 end
% %             end
% %         end
% %         
% %         fprintf('Calculating SLOs...\n');
% %         for i_channel = 1:size(LFP_mV_adj,2)
% %             if sum(LFP_mV_adj{1,i_channel}) == 0
% %                 SLOcount(:,i_channel) = NaN;
% %             else
% %                 for i_trial = 1:size(LFP_mV_adj,1)
% %                     mV = LFP_mV_adj{i_trial,i_channel};
% %                     
% %                     
% %                     % Filter the data
% %                     mV = filtfilt(b, a, mV); % BP filter
% %                     
% %                     mV = abs(mV);
% %                     thresh = mean(mV) + 6*std(mV);
% %                     over_idx = find(mV > thresh);
% %                     over_idx = over_idx(over_idx > 20000 & over_idx <= 40000); % Only look during stim
% %                     numSLO = 0;
% %                     idxSLO = [];
% %                     for i_over = 1:length(over_idx)
% %                         
% %                         % calculate indices for 500 ms before/after current time point
% %                         start = max(1,over_idx(i_over)-minT);
% %                         stop = min(over_idx(i_over)+postT, length(mV));
% %                         
% %                         % extract time series segments corresponding to 500 ms before/after
% %                         % current time point
% %                         pre = mV(start:(over_idx(i_over)-1));
% %                         post = mV(over_idx(i_over):stop);
% %                         
% %                         % does current time point satisfy the algorithm requirements? if so,
% %                         % count it as the beginning of a SLO.
% %                         if all(pre < thresh) && mean(post>thresh)>=postThresh
% %                             numSLO = numSLO+1;
% %                             idxSLO = [idxSLO over_idx(i_over)];
% %                         end
% %                     end
% %                     SLOcount(i_trial,i_channel) = numSLO;
% %                     SLOidx{i_trial,i_channel} = idxSLO;
% %                 end
% %             end
% %         end
% %         
% %         % Save the results
% %         fprintf('Saving: %s  BF_%s_%sHz_countSLOs.mat\n\n',mouse,region,frequency);
% %         SaveFilename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_countSLOs.mat',region,frequency));
% %         save(SaveFilename,'SLOcount','SLOidx');
% %         
% %     end
% % end

%%

minT = 500; % Rest period before SLO (ms)
postT = 500; % Duriation of SLO measurement (ms)
postThresh = 0.02; % Threshold for SLO during postT

fc1 = 2; % highpass frequency
fc2 = 200; % lowpass frequency
fs = 1000;
Wn = [fc1 fc2]*2/fs;
ftype = 'bandpass';
[b,a] = butter(5, Wn, ftype);

for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
    if exist(Filename,'file')
        fprintf('%d - Loading: %s\n',i,Filename);
        load(Filename);
        
        SLOcount_pre = zeros(size(LFP_mV_adj));
        SLOcount_peri = zeros(size(LFP_mV_adj));
        SLOcount_post = zeros(size(LFP_mV_adj));
        SLOidx_pre = cell(size(LFP_mV_adj));
        SLOidx_peri = cell(size(LFP_mV_adj));
        SLOidx_post = cell(size(LFP_mV_adj));
        
                
        for xx = 1:size(LFP_mV_adj,1)
            for yy = 1:size(LFP_mV_adj,2)
                if isempty(LFP_mV_adj{xx,yy})
                    LFP_mV_adj{xx,yy} = zeros(1,60000);
                end
            end
        end
        
        fprintf('Calculating SLOs...\n');
        for i_channel = 1:size(LFP_mV_adj,2)
            if sum(LFP_mV_adj{1,i_channel}) == 0
                SLOcount_pre(:,i_channel) = NaN;
                SLOcount_peri(:,i_channel) = NaN;
                SLOcount_post(:,i_channel) = NaN;
            else
                for i_trial = 1:size(LFP_mV_adj,1)
                    mV = LFP_mV_adj{i_trial,i_channel};
                    
                    
                    % Filter the data
                    mV = filtfilt(b, a, mV); % BP filter
                    
                    mV = abs(mV);
                    thresh = mean(mV) + 6*std(mV);
                    over_idx_full = find(mV > thresh);
                    
                    
                    over_idx = over_idx_full(over_idx_full <= 20000); % Pre-stimulation
                    numSLO = 0;
                    idxSLO = [];
                    for i_over = 1:length(over_idx)
                        
                        % calculate indices for 500 ms before/after current time point
                        start = max(1,over_idx(i_over)-minT);
                        stop = min(over_idx(i_over)+postT, length(mV));
                        
                        % extract time series segments corresponding to 500 ms before/after
                        % current time point
                        pre = mV(start:(over_idx(i_over)-1));
                        post = mV(over_idx(i_over):stop);
                        
                        % does current time point satisfy the algorithm requirements? if so,
                        % count it as the beginning of a SLO.
                        if sum(pre < thresh) == minT && sum(post>thresh) >= (postThresh*postT)
                            numSLO = numSLO+1;
                            idxSLO = [idxSLO over_idx(i_over)];
                        end
                    end
                    SLOcount_pre(i_trial,i_channel) = numSLO;
                    SLOidx_pre{i_trial,i_channel} = idxSLO;
                    
                    over_idx = over_idx_full(over_idx_full > 20000 & over_idx_full <= 40000); % Peri-stimulation
                    numSLO = 0;
                    idxSLO = [];
                    for i_over = 1:length(over_idx)
                        
                        % calculate indices for 500 ms before/after current time point
                        start = max(1,over_idx(i_over)-minT);
                        stop = min(over_idx(i_over)+postT, length(mV));
                        
                        % extract time series segments corresponding to 500 ms before/after
                        % current time point
                        pre = mV(start:(over_idx(i_over)-1));
                        post = mV(over_idx(i_over):stop);
                        
                        % does current time point satisfy the algorithm requirements? if so,
                        % count it as the beginning of a SLO.
                        if sum(pre < thresh) == minT && sum(post>thresh) >= (postThresh*postT)
                            numSLO = numSLO+1;
                            idxSLO = [idxSLO over_idx(i_over)];
                        end
                    end
                    SLOcount_peri(i_trial,i_channel) = numSLO;
                    SLOidx_peri{i_trial,i_channel} = idxSLO;
                    
                    over_idx = over_idx_full(over_idx_full > 40000 & over_idx_full <= 60000); % Post-stimulation
                    numSLO = 0;
                    idxSLO = [];
                    for i_over = 1:length(over_idx)
                        
                        % calculate indices for 500 ms before/after current time point
                        start = max(1,over_idx(i_over)-minT);
                        stop = min(over_idx(i_over)+postT, length(mV));
                        
                        % extract time series segments corresponding to 500 ms before/after
                        % current time point
                        pre = mV(start:(over_idx(i_over)-1));
                        post = mV(over_idx(i_over):stop);
                        
                        % does current time point satisfy the algorithm requirements? if so,
                        % count it as the beginning of a SLO.
                        if sum(pre < thresh) == minT && sum(post>thresh) >= (postThresh*postT)
                            numSLO = numSLO+1;
                            idxSLO = [idxSLO over_idx(i_over)];
                        end
                    end
                    SLOcount_post(i_trial,i_channel) = numSLO;
                    SLOidx_post{i_trial,i_channel} = idxSLO;
                end
            end
        end
        
        % Save the results
        fprintf('Saving: %s  BF_%s_%sHz_countSLOs_PrePeriPost.mat\n\n',mouse,region,frequency);
        SaveFilename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_countSLOs_PrePeriPost.mat',region,frequency));
        save(SaveFilename,'SLOcount_pre','SLOidx_pre','SLOcount_peri','SLOidx_peri','SLOcount_post','SLOidx_post');
        
        SLOcount = SLOcount_peri;
        SLOidx = SLOidx_peri;
        fprintf('Saving: %s  BF_%s_%sHz_countSLOs.mat\n\n',mouse,region,frequency);
        SaveFilename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_countSLOs.mat',region,frequency));
        save(SaveFilename,'SLOcount','SLOidx');
    end
end