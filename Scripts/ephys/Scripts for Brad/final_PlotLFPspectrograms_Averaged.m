% Plot and save the LFP Spectrogram averaged across trials

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

fileTag = '_filt';

t_window = 1000;
t_overlap = 800;
nfft = 500;
fs = 1000;
climit = [-75 -35];
ylimit = [0 80];

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
mouseList(:,5) = cellfun(@(x) double(strsplit(string(x),',')),raw(rows,21),'UniformOutput',false);
mouseList(:,6) = cellfun(@(x) double(strsplit(string(x),',')),raw(rows,22),'UniformOutput',false);
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        
        n_channels = size(LFP_mV_adj,2);
        n_trials = size(LFP_mV_adj,1);
        
        % Calculate the averages
        fprintf('Averaging LFP...\n');
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
                end
                trials_per_channel = reshape(cell2mat(LFP_mV_adj(:,i_channel)),length(LFP_mV_adj{1,1}),n_trials)';
                average_trials{1,i_channel} = mean(trials_per_channel,1);
            else
                average_trials{1,i_channel} = zeros(1,60000);
            end
        end
        
        % Plot the spectrograms of the averages
        fprintf('Plotting/Saving...\n');
        for i_channel = 1:n_channels
            p = 16;
            subplot(4,4,mod(i_channel-1,p)+1);
            spectrogram(average_trials{1,i_channel},t_window,t_overlap,nfft,fs,'yaxis');
            hold on;
            text(0.02,0.90,['Ch.' char(num2str(i_channel))],'Units','Normalized','FontWeight','bold','Color','w');
            % Format the plots
            ax=gca;
            ax.CLim = climit;
            ylim(ylimit);
            
            % Save the plot
            if mod(i_channel-1,p)+1 == p
                xlabel(sprintf('%s       %s       %sHz       %s       Ch.%d-%d',group,region,frequency,mouse,i_channel-15,i_channel),'Interpreter','none');
                saveas(Fig1,fullfile(working_dir,'Figures_Indiv','PlotLFPspectrograms_Averaged',sprintf('%s_%s_%sHz_E%d.png',mouse,region,frequency,i_channel/16)),'png');
%                 pause;
                clf;
            end
        end
    end
end