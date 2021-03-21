function final_CheckLFPtraces_BE(storage,base_fold,slash)

% Plot the LFP traces (modified to save individual plots)

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

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
stim = {'0_1' '0_5' '1_0'};
for i_mouse = 1:size(base_fold,1)
    
    SLOcount = zeros(size(stim,2),32);  
    SLOidx = cell(size(stim,2),32);
    
    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '.mat'];
        save_file = [data_fold stim{i_stim} '_countSLOS.mat'];
        
        if exist(data_file,'file')
            load(data_file);
                    
            TS_delta = Laser_TS(1) - 24;
            mincut = 0; maxcut = 600;
            if TS_delta < 0
                Diff = find(diff(Laser_TS) > 1);
                Diff = Diff(1) + 1;
                TS_delta = Laser_TS(Diff) - 24;
                mincut = 60; maxcut = 660;
            end
            cutoffs = [mincut:60:maxcut] + TS_delta;
            offset = cutoffs(1)*1000;
            idx_offset = offset;
                    
            % Plot the LFP traces
            for i_ch = 1:32
                p = 32;
                figure(Fig1);subplot(8,4,mod(i_ch-1,p)+1);

                % Find possible clipping
                mV = LFP_mV{i_ch};
                max_mV = max(mV); min_mV = min(mV); lim_mV = max(max_mV,abs(min_mV));
                max_mV = lim_mV; min_mV = -lim_mV;
                max_idx = find(mV == max_mV);
                min_idx = find(mV == min_mV);
                if (length(max_idx)+length(min_idx)) == 1
                    max_idx = []; min_idx = [];
                end
                max_pts = LFP_TS(max_idx);
                min_pts = LFP_TS(min_idx);

                % Filter the data
                filt_mV = filtfilt(b, a, mV); % BP filter
                % Find the SLOs
                temp_TS = abs(filt_mV);
                thresh = mean(temp_TS) + 6*std(temp_TS);
                
                over_idx = find(temp_TS > thresh);
                over_idx2 = [];
                for k = 1:10
                    over_idx2 = [over_idx2; over_idx(over_idx > (20000 + 60000*(k-1) + offset) & over_idx < (40000 + 60000*(k-1) + offset))];
                end
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

                SLOcount(i_stim,i_ch) = numSLO;
                SLOidx{i_stim,i_ch} = idxSLO;

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
                set(gca,'Position', pos);
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
                xlim([(cutoffs(1)-5) cutoffs(end)+5]);
                ylim([min_mV-0.3 max_mV+0.3]);
            end
            % Save the plot
            saveas(Fig1,[data_fold 'CheckLFPtraces_' stim{i_stim} '.png']);
            
            fprintf('\n')
            fprintf('Saving: %s',  save_file);
            fprintf('\n')
            save(save_file,'SLOcount','SLOidx','idx_offset');
            clf
        end
    end
end
close all;