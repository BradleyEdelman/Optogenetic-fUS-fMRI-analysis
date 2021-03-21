% Plot and save the LFP Spectrogram averaged across mice

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
namelist = {'low','high'};

fileTag = '_filt';

t_window = 1000;
t_overlap = 800;
nfft = 500;
fs = 1000;
climit = [-75 -35];
ylimit = [0 80];

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');

fprintf('Loading: LFPaverages.mat\n');
Filename = fullfile(working_dir,'LFPaverages.mat');

load(Filename);

for i_group = 1:length(Groups)
    group = Groups{i_group};
    for i_region = 1:length(Regions)
        region = Regions{i_region};
        if i_region == 1
            d = 0;
        else
            d = 16;
        end
        for i_frequency = 1:size(All_Data,3)
            if ~isempty(All_Data{i_group,i_region,1,1})
                
                fprintf('Plotting/Saving: %s %s %s\n',group,region,namelist{i_frequency});
                for i_channel = 1:16
                    subplot(4,4,i_channel);
                    temp_data = mean(All_Data{i_group,i_region,i_frequency,i_channel},1);
                    spectrogram(temp_data,t_window,t_overlap,nfft,fs,'yaxis');
                    hold on;
                    text(0.02,0.90,['Ch.' char(num2str(i_channel+d))],'Units','Normalized','FontWeight','bold','Color','w');
                    % Format the plots
                    ax=gca;
                    ax.CLim = climit;
                    ylim(ylimit);
                end

                % Save the plot
                xlabel(sprintf('%s       %s       %s',group,region,namelist{i_frequency}),'Interpreter','none');
                saveas(Fig1,fullfile(working_dir,'Figures_Group','PlotLFPspectrograms_byGroup',sprintf('LFPspectrogram_%s_%s_%s.png',region,group,namelist{i_frequency})),'png');
%                 pause;
                clf;
            end
        end
    end
end