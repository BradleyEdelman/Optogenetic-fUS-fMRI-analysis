function final_PlotLFP_Info_BE(storage,base_fold,slash,type)
% Plot LFP spectrogram and amplitude change across animals

t_window = 1000;
t_overlap = 800;
nfft = 500;
fs = 1000;
climit = [-75 -20];
ylimit = [0 80];

%%
Fsp(1) = figure('Position', get(0, 'Screensize'),'visible','off');
Fsp(2) = figure('Position', get(0, 'Screensize'),'visible','off');
Fsp(3) = figure('Position', get(0, 'Screensize'),'visible','off');
Fsp_region(1) = figure('Position', get(0, 'Screensize'),'visible','off');
Fsp_region(2) = figure('Position', get(0, 'Screensize'),'visible','off');
Fsp_region(3) = figure('Position', get(0, 'Screensize'),'visible','off');

Ftot(1) = figure('Position', get(0, 'Screensize'),'visible','off');
Ftot(2) = figure('Position', get(0, 'Screensize'),'visible','off');

FLFP(1) = figure('Position', get(0, 'Screensize'),'visible','off');
FLFP(2) = figure('Position', get(0, 'Screensize'),'visible','off');

Famp = figure('Position', get(0, 'Screensize'),'visible','off');

stim = {'0_1' '0_5' '1_0'};

for i_stim = 1:size(stim,2)
    stim_fold = [storage stim{i_stim} slash];
    if ~exist(stim_fold,'dir'); mkdir(stim_fold); end
    
    for i_mouse = 1:size(base_fold,1)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        clean_file = [data_fold stim{i_stim} '_clean.mat'];
        LFP_file = [data_fold stim{i_stim} '_LFP_amp.mat'];
%         save_file = [data_fold stim{i_stim} '_LFP.mat'];
        
        if exist(clean_file,'file') && exist(LFP_file,'file')
            
            load(clean_file)
            
            n_channels = size(LFP_mV_adj,2);
            n_trials = size(LFP_mV_adj,1);
        
            % Store each trial per animal
            for i_channel = 1:n_channels
                if ~isempty(vertcat(LFP_mV_adj{:,i_channel}))
                    for i_trial = 1:n_trials
                        if isempty(LFP_mV_adj{i_trial,i_channel})
                        elseif length(LFP_mV_adj{i_trial,i_channel}) > 60000
                            LFP_mV_adj{i_trial,i_channel} = LFP_mV_adj{i_trial,i_channel}(1:60000);
                        elseif length(LFP_mV_adj{i_trial,i_channel}) < 60000
                            endpoint = length(LFP_mV_adj{i_trial,i_channel});
                            LFP_mV_adj{i_trial,i_channel}((endpoint+1):60000) = LFP_mV_adj{i_trial,i_channel}(endpoint);
                        end
                    end
                    trials_per_channel = horzcat(LFP_mV_adj{:,i_channel})';
                    
                    if n_trials > 10
                        trials_per_channel(1,:) = [];
                    end
                    
                    if ~isempty(trials_per_channel)
                        average_trials{i_mouse,i_channel} = nanmean(trials_per_channel,1);
                    else
                        average_trials{i_mouse,i_channel} = nan(1,60000);
                    end
                else
                    average_trials{i_mouse,i_channel} = nan(1,60000);
                end
            end
            
            average_trials_region{i_mouse,1,i_stim} = nanmean(vertcat(average_trials{i_mouse,1:8}),1);
            average_trials_region{i_mouse,2,i_stim} = nanmean(vertcat(average_trials{i_mouse,9:16}),1);
            average_trials_region{i_mouse,3,i_stim} = nanmean(vertcat(average_trials{i_mouse,17:24}),1);
            average_trials_region{i_mouse,4,i_stim} = nanmean(vertcat(average_trials{i_mouse,25:32}),1);
            
            % Store each region per animal
            load(LFP_file)
            LFP_regions_amp(i_mouse,:) =...
                [nanmean(LFP_amp_region{1}(:,2)) nanmean(LFP_amp_region{2}(:,2)) nanmean(LFP_amp_region{3}(:,2)) nanmean(LFP_amp_region{4}(:,2))];
            LFP_regions_amp_change(i_mouse,:) = LFP_amp_region_change(1,:);
            LFP_z_regions(i_mouse,:) = LFP_ave_z_region(1,:);
            
        end
        
    end

    % Average each LFP trace for each channel, plot average spectrogram (Group)
    for i_channel = 1:size(average_trials,2)
        animal_trial_ts_GRP(i_channel,:,i_stim) = mean(vertcat(average_trials{:,i_channel}),1);
        
%         figure(Fsp(i_stim)); subplot(4,8,i_channel);
%         spectrogram(animal_trial_ts_GRP(i_channel,:,i_stim),t_window,t_overlap,nfft,fs,'yaxis');
%         hold on;
%         text(0.02,0.90,['Ch.' char(num2str(i_channel))],'Units','Normalized','FontWeight','bold','Color','w');
%         % Format the plots
%         ax=gca; ax.CLim = climit;
%         ylim(ylimit);
    end
    
    % Plot the spectrograms of the averages per region (Group)
    average_trial_ts_region_GRP{1,i_stim} = nanmean(animal_trial_ts_GRP(1:8,:,i_stim),1);
    average_trial_ts_region_GRP{2,i_stim} = nanmean(animal_trial_ts_GRP(9:16,:,i_stim),1);
    average_trial_ts_region_GRP{3,i_stim} = nanmean(animal_trial_ts_GRP(17:24,:,i_stim),1);
    average_trial_ts_region_GRP{4,i_stim} = nanmean(animal_trial_ts_GRP(25:32,:,i_stim),1);
    t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
%     for i_region = 1:4
%         figure(Fsp_region(i_stim)); subplot(2,2,i_region);
%         spectrogram(average_trial_ts_region_GRP{i_region,i_stim},t_window,t_overlap,nfft,fs,'yaxis');
%         hold on;
%         text(0.02,0.90,[t_bank{i_region}],'Units','Normalized','FontWeight','bold','Color','w');
%         % Format the plots
%         ax=gca; ax.CLim = climit;
%         ylim(ylimit);
%     end
    
    LFP_region_GRP_mV{i_stim} = LFP_regions_amp;
    LFP_region_GRP{i_stim} = LFP_regions_amp_change;
    LFP_z_region_GRP{i_stim} = LFP_z_regions;   
end


% Plot individual animal LFP traces per region and stim intensity
t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
for i_mouse = 1:size(base_fold,1)
    C = {'b','r','g'};
    figure(FLFP(1)); clf; figure(FLFP(2)); clf
    for i_region = 1:4
        for i_stim = size(stim,2):-1:1
            figure(FLFP(1)); subplot(3,4,i_region + 4*(i_stim-1)); hold on
            
%             fc1 = 8; % highpass frequency
%             fc2 = 12; % lowpass frequency
%             Wn = [fc1 fc2]*2/fs;
%             [b,a] = butter(4, Wn, 'bandpass');
%             
%             tmp = filtfilt(b,a,average_trials_region{i_mouse,i_region,i_stim});
%             figure; plot(tmp)
%             
%             
%             
            
            [pks,locs] = findpeaks(average_trials_region{i_mouse,i_region,i_stim},1000);
            locs(pks<0) = []; pks(pks<0) = []; 
            [pks2,locs2] = findpeaks(-average_trials_region{i_mouse,i_region,i_stim},1000); pks2 = -pks2;
            locs2(pks2>0) = []; pks2(pks2>0) = []; 
            
            [y, ty] = resample([0 pks 60],[0 locs 60],1000); y = y(1500:58500);
            [y2, ty2] = resample([0 pks2 60],[0 locs2 60],1000); y2 = y2(1500:58500);
            
            plot(y); plot(y2)
            ylim([-4.5 2])
            title([t_bank{i_region} ':' stim{i_stim}])
            
            figure(FLFP(2)); subplot(2,2,i_region); hold on
            plot(average_trials_region{i_mouse,i_region,i_stim});
            ylim([-4.5 2])
            title([t_bank{i_region}])
        end
    end
    figure(FLFP(1)); supertitle(['Mouse: ' num2str(i_mouse)])
    figure(FLFP(2)); supertitle(['Mouse: ' num2str(i_mouse)])
    
    data_fold = [storage base_fold{i_mouse} slash];
    saveas(FLFP(1),[data_fold 'LFP_envelope_total.svg']);
    saveas(FLFP(2),[data_fold 'LFP_trace_total.svg']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% GROUP RESULTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% 
% Average LFP Trace across animals
for i_region = 1:4
    for i_stim = 3:-1:1
    
        figure(Ftot(1)); subplot(2,2,i_region); hold on
        plot(average_trial_ts_region_GRP{i_region,i_stim});
        
        % Format the plots
        ylim([-1.5 1.5]);
        ylim([-2.5 2])
        title([t_bank{i_region}])
    end
end
supertitle('Mouse: Group Average')
saveas(Ftot(1),[storage type '_LFP_envelope_total_GROUP.svg']);


% Average spectrogram across animals
for i_stim = 1:3    
    for i_region = 1:4
    
        figure(Ftot(2)); subplot(3,4,i_region + 4*(i_stim-1));
        spectrogram(average_trial_ts_region_GRP{i_region,i_stim},t_window,t_overlap,nfft,fs,'yaxis');
        hold on;
        text(0.02,0.90,[stim{i_stim} ':' t_bank{i_region}],'Units','Normalized','FontWeight','bold','Color','w');
        % Format the plots
        ax=gca; ax.CLim = climit;
        ylim(ylimit);
        
        LFPAVE(i_stim,i_region) = nanmean(LFP_region_GRP{i_stim}(:,i_region),1);
        LFPsem(i_stim,i_region) = nanstd(LFP_region_GRP{i_stim}(:,i_region),1)/sqrt(size(LFP_region_GRP{i_stim},1));
    end
end
supertitle('Mouse: Group Average')
saveas(Ftot(1),[storage type '_Spectrogram_total_GROUP.svg']);
save([storage type '_LFP_AMP_values_GROUP'],'LFP_region_GRP','LFP_region_GRP_mV','LFP_z_region_GRP','base_fold');

% bar plot percent change
figure; hold on
ctrs = 1:4; hBar{1} = bar(ctrs, LFPAVE'); drawnow; clear ctr ydt
for k1 = 1:size(LFPAVE',2)
    ctr(k1,:) = bsxfun(@plus, hBar{1}(1).XData, [hBar{1}(k1).XOffset]');
    ydt(k1,:) = hBar{1}(k1).YData;
end
errorbar(ctr, ydt, zeros(4,3)', LFPsem, zeros(4,3)', zeros(4,3)', '.k',...
    'capsize',5,'marker','o','markersize',0.1,'linewidth',1)
hBar{2} = bar(ctrs, LFPAVE'); delete(hBar{1})
set(gca,'ylim',[0 1e4],'xtick',[1 2 3 4],'xticklabel',t_bank);
saveas(Famp,[storage type '_LPF_AMP_total_GROUP.svg']);

% box plot percent change
figure;
for k2 = 1:4
    subplot(2,2,k2)
    D = [LFP_region_GRP{1}(:,k2);...
        LFP_region_GRP{2}(:,k2);...
        LFP_region_GRP{3}(:,k2)];
    D1 = [repmat('1',4,1); repmat('2',4,1); repmat('3',4,1)];
    boxplot(D, D1)
    
    title(t_bank{k2})
    set(gca,'ylim',[0 1.5e4]);
end




            