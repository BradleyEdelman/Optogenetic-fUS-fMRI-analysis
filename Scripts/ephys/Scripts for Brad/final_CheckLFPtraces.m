% Plot the LFP traces (modified to save individual plots)

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';
Mice = {'061519_A','061519_B','061519_C','061619_C','061719_A','061719_B','061719_C','061819_A','061819_B','061819_C','061819_D','062019_A','062019_B','062019_C','062019_D','062119_A','062119_B','063019_A','063019_B','063019_C','070119_A','070119_B','070119_C','070619_A','070619_B','070619_C','070719_A','070719_B','070719_C','070719_D','070819_A','070819_B','070819_C'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'6','10','20'};
Recordings = {'A','B','C','D'};
TrialStart = {1,6,11,16};
TrialEnd = {5,10,15,20};
i_png = 529;

minT = 500; % Rest period before SLO (ms)
postT = 500; % Duriation of SLO measurement (ms)
postThresh = 0.02; % Threshold for SLO during postT

fc1 = 2; % highpass frequency
fc2 = 200; % lowpass frequency
fs = 1000;
Wn = [fc1 fc2]*2/fs;
ftype = 'bandpass';
[b,a] = butter(5, Wn, ftype);

Mice = {'072419_A'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
for i_mouse = 1:length(Mice)
    mouse = Mice{i_mouse};
    for i_region = 1:length(Regions)
        region = sprintf('BF_%s',Regions{i_region});
        for i_freq = 1:length(Frequencies)
            frequency = Frequencies{i_freq};
            SLOcount = zeros(4,32);
            SLOidx = cell(4,32);
            save_flag = 0;
            for i_recording = 1:length(Recordings)
                recording = Recordings{i_recording};
                Filename = fullfile(working_dir,mouse,sprintf('%s_%sHz_%s.mat',region,frequency,recording));
                if exist(Filename,'file')
                    fprintf('Loading: %s\n',Filename);
                    load(Filename);
                    save_flag = 1;
                    
                    TS_delta = Laser_TS(1) - 20;
                    cutoffs = [0:60:300] + TS_delta;
                    offset = cutoffs(1)*1000;
                    idx_offset{i_recording} = offset;
                    
                    % Plot the LFP traces
                    fprintf('Plotting/Saving: LFP Traces...\n');
                    for i_ch = 1:32
                        p = 32;
                        subplot(8,4,mod(i_ch-1,p)+1);
                        
                        % Find possible clipping
                        mV = LFP_mV{i_ch};
                        max_mV = max(mV);
                        min_mV = min(mV);
                        lim_mV = max(max_mV,abs(min_mV));
                        max_mV = lim_mV;
                        min_mV = -lim_mV;
                        max_idx = find(mV == max_mV);
                        min_idx = find(mV == min_mV);
                        if (length(max_idx)+length(min_idx)) == 1
                            max_idx = [];
                            min_idx = [];
                        end
                        max_pts = LFP_TS(max_idx);
                        min_pts = LFP_TS(min_idx);
                        
                        % Filter the data
                        filt_mV = filtfilt(b, a, mV); % BP filter
                        % Find the SLOs
                        temp_TS = abs(filt_mV);
                        thresh = mean(temp_TS) + 6*std(temp_TS);
                        over_idx = find(temp_TS > thresh);
                        over_idx_1 = over_idx(over_idx > (20000+offset) & over_idx < (40000+offset));
                        over_idx_2 = over_idx(over_idx > (80000+offset) & over_idx < (100000+offset));
                        over_idx_3 = over_idx(over_idx > (140000+offset) & over_idx < (160000+offset));
                        over_idx_4 = over_idx(over_idx > (200000+offset) & over_idx < (200000+offset));
                        over_idx_5 = over_idx(over_idx > (260000+offset) & over_idx < (260000+offset));
                        over_idx = [over_idx_1' over_idx_2' over_idx_3' over_idx_4' over_idx_5'];
                        numSLO = 0;
                        idxSLO = [];
                        for i_over = 1:length(over_idx)
                            % calculate indices for 500 ms before/after current time point
                            start = max(1,over_idx(i_over)-minT);
                            stop = min(over_idx(i_over)+postT, length(temp_TS));
                            
                            % extract time series segments corresponding to 500 ms before/after
                            % current time point
                            pre = temp_TS(start:(over_idx(i_over)-1));
                            post = temp_TS(over_idx(i_over):stop);
                            
                            % does current time point satisfy the algorithm requirements? if so,
                            % count it as the beginning of a SLO.
                            if all(pre < thresh) && mean(post>thresh)>postThresh
                                numSLO = numSLO+1;
                                idxSLO = [idxSLO over_idx(i_over)];
                            end
                        end
                        
                        SLOcount(i_recording,i_ch) = numSLO;
                        SLOidx{i_recording,i_ch} = idxSLO;
                        
                        % Plot the results
                        hold on;
                        for i_cutoff = 1:length(cutoffs)
                            line([cutoffs(i_cutoff) cutoffs(i_cutoff)],[min_mV-0.3,max_mV+0.3],'Color','m');
                        end
                        plot1 = plot(LFP_TS,mV,'r');
                        plot1.Color(4) = 0.3;
                        plot(LFP_TS,filt_mV,'k');
                        scatter(Laser_TS,(min_mV-0.2)*ones(length(Laser_TS),1),'.b');
                        scatter(LFP_TS(idxSLO),(max_mV+0.1)*ones(length(idxSLO),1),'filled','MarkerFaceColor',[0 0.6700 0.1900]);
                        scatter(max_pts,(max_mV+0.1)*ones(length(max_pts),1),'vr');
                        scatter(min_pts,(min_mV-0.1)*ones(length(min_pts),1),'^r');
                        text(cutoffs(1)+5,max_mV+0.15,[char(num2str(i_ch))]);
%                         text(cutoffs(end)-30,max_mV+0.15,['SLOs.' char(num2str(length(SLO_markers{6,i_ch})))])
%                         text(cutoffs(1)+5,min_mV+-.10,['Clips.' char(num2str(length(max_idx)+length(min_idx)))])
                        
                        % Format the plots
                        pos = get(gca, 'Position');
                        pos(1) = 0.003+0.250*(mod(i_ch-1,4));
                        pos(3) = 0.245;
                        pos(4) = 0.1;
                        set(gca, 'Position', pos);
                        set(gca,'XTick',[]);
                        set(gca,'YTick',[]);
                        xlim([(cutoffs(1)-5) cutoffs(end)+5]);
                        ylim([min_mV-0.3 max_mV+0.3]);
                        
                        % Save the plot
                        if mod(i_ch-1,p)+1 == p
                            xlabel(sprintf('%s %s %sHz Trials %d-%d',mouse,region(4:end),frequency,TrialStart{i_recording},TrialEnd{i_recording}),'Interpreter','none');
                            saveas(Fig1,fullfile(working_dir,'Figures_Indiv','CheckLFPtraces',sprintf('%d.png',i_png)),'png');
                            i_png = i_png+1;
%                             pause;
                            clf;
                        end
                    end
                end
            end
            if save_flag == 1
                % Save the SLOs
                fprintf('Saving: %s  %s_%sHz_countSLOs_prelim.mat\n\n',mouse,region,frequency);
                SaveFilename = fullfile(working_dir,mouse,sprintf('%s_%sHz_countSLOs.mat',region,frequency));
                save(SaveFilename,'SLOcount','SLOidx','idx_offset');
            end
        end
    end
end
close all;