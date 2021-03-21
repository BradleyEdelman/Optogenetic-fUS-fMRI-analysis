% Plot the SLOs

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

TrialStart = {1,6,11,16};
TrialEnd = {5,10,15,20};
i_png = 1;

fc1 = 2; % highpass frequency
fc2 = 200; % lowpass frequency
fs = 1000;
Wn = [fc1 fc2]*2/fs;
ftype = 'bandpass';
[b,a] = butter(5, Wn, ftype);

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','on');
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_filt.mat',region,frequency));
    Filename2 = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_countSLOs.mat',region,frequency));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        load(Filename2);
        
        % Fill in missing channels & trials
        if ~isempty(badChannels)
            for i_BC = badChannels
                for i_trial = 1:size(LFP_mV_adj,1)
                    LFP_mV_adj{i_trial,i_BC} = zeros(60000,1);
                    SLOcount(i_trial,i_BC) = 0;
                    SLOidx{i_trial,i_BC} = [];
                end
            end
        end
        
        if ~isempty(badTrials)
            for i_ch = 1:32
                bump = length(badTrials);
                for i_BT = 20:-1:1
                    if ~ismember(i_BT,badTrials)
                        Laser_TS_adj{1,i_BT} = Laser_TS_adj{1,i_BT - bump};
                        LFP_TS_adj{1,i_BT} = LFP_TS_adj{1,i_BT - bump};
                        LFP_mV_adj{i_BT,i_ch} = LFP_mV_adj{i_BT - bump,i_ch};
                        SLOcount(i_BT,i_ch) = SLOcount(i_BT - bump,i_ch);
                        SLOidx{i_BT,i_ch} = SLOidx{i_BT - bump,i_ch};
                    else
                        Laser_TS_adj{1,i_BT} =  Laser_TS_adj{1,1};
                        LFP_TS_adj{1,i_BT} = LFP_TS_adj{1,1};
                        LFP_mV_adj{i_BT,i_ch} = zeros(60000,1);
                        SLOcount(i_BT,i_ch) = 0;
                        SLOidx{i_BT,i_ch} = [];
                        bump = bump - 1;
                    end
                end
            end
        end       
        
        % Plot the traces with the SLOs
        fprintf('Plotting/Saving: LFP Traces...\n');
        
        for i_trialblock = 1:length(TrialStart)
            
            for i_ch = 1:32
                p = 32;
                subplot(8,4,mod(i_ch-1,p)+1);
                
                for i_trial = TrialStart{i_trialblock}:TrialEnd{i_trialblock}
                    delta = 60*mod(i_trial-1,5);
                    
                    LFP_TS = LFP_TS_adj{1,i_trial} + delta;
                    Laser_TS = Laser_TS_adj{1,i_trial} + delta;
                    idxSLO = SLOidx{i_trial,i_ch};
                    LFP_mV = LFP_mV_adj{i_trial,i_ch};
                    
                    if length(LFP_mV) > 60000
                        LFP_mV = LFP_mV(1:60000);
                    end
                    if length(LFP_mV) < length(LFP_TS)
                        LFP_TS = LFP_TS(1:length(LFP_mV));
                    end

                    % Find the limits
                    if delta == 0
                        lim_mV = max([abs(LFP_mV_adj{i_trial,i_ch});abs(LFP_mV_adj{i_trial+1,i_ch});abs(LFP_mV_adj{i_trial+2,i_ch});abs(LFP_mV_adj{i_trial+3,i_ch});abs(LFP_mV_adj{i_trial+4,i_ch})]);
                        max_mV = lim_mV;
                        min_mV = -lim_mV;
                    end
         
                    % Filter the data
                    filt_mV = filtfilt(b, a, LFP_mV); % BP filter
                    
                    % Plot the data
                    hold on;
                    line([delta delta],[min_mV-0.3,max_mV+0.3],'Color','m');
                    plot1 = plot(LFP_TS,LFP_mV,'r');
                    plot1.Color(4) = 0.3;
                    plot(LFP_TS,filt_mV,'k');
                    scatter(Laser_TS,(min_mV-0.2)*ones(length(Laser_TS),1),'.b');
                    scatter(LFP_TS(idxSLO),(max_mV+0.1)*ones(length(idxSLO),1),'filled','MarkerFaceColor',[0 0.6700 0.1900]);
                end
                text(5,max_mV+0.15,[char(num2str(i_ch))]);
                
                % Format the plots
                pos = get(gca, 'Position');
                pos(1) = 0.003+0.250*(mod(i_ch-1,4));
                pos(3) = 0.245;
                pos(4) = 0.1;
                set(gca, 'Position', pos);
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
                xlim([-5 305]);
                ylim([min_mV-0.3 max_mV+0.3]);
                
            end
            % Save the plot
            if mod(i_ch-1,p)+1 == p
                xlabel(sprintf('%s %s %sHz Trials %d-%d',mouse,region,frequency,TrialStart{i_trialblock},TrialEnd{i_trialblock}),'Interpreter','none');
                fprintf('Saving: %s\n',fullfile(working_dir,'Figures_Indiv','PlotSLOtraces',sprintf('%d.png',i_png)));
%                 saveas(Fig1,fullfile(working_dir,'Figures_Indiv','PlotSLOtraces',sprintf('%d.png',i_png)),'png');
                pause;
                clf;
                i_png = i_png+1;
            end
        end
    end
end