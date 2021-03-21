% Plot the firing rates per period per neuron for all neurons colored by
% response type by electrode site

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

hist_binsize = .00025; %seconds
xlimit = [0 .03];

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_SpikeDelays.mat',region,frequency));
    Filename2 = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_spikeTest.mat',region,frequency));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        load(Filename2);
        
        num_trials = size(Spike_Delays,1);
        num_channels = size(Spike_Delays,2);
        num_units = size(Spike_Delays,3);
                
        % Loop through and plot all of the neurons' Spike Delays
        for i_ch = 1:num_channels
            p = 16;
            subplot(4,4,mod(i_ch-1,p)+1);
            hold on;
            
            for i_unit = 1:num_units
                if ~all(cellfun('isempty',Spike_Delays(:,i_ch,i_unit)))
                    if INC(1,i_ch,i_unit) == 1
                        barcolor = [0.6350 0.0780 0.1840];
                    elseif DEC(1,i_ch,i_unit) == 1
                        barcolor = [0 0.4470 0.7410];
                    elseif NC(1,i_ch,i_unit) == 1
                        barcolor = [0.4660 0.6740 0.1880];
                    end
                    pts = cell2mat(Spike_Delays(:,i_ch,i_unit));
                    pts = pts(pts >= 1/40000); % Get rid of values less than the sampling frequency
                    histogram(pts,[0:hist_binsize:0.200],'FaceColor',barcolor,'FaceAlpha',0.4);
                    hold on;
                end
            end
            text(0.02,0.90,['Ch.' char(num2str(i_ch))],'Units','Normalized');
            % Format the plots
            xticks = [0:0.005:0.20];
            set(gca,'XTick',xticks);
            set(gca,'XTickLabel',xticks*1000);
            set(gca,'TickDir','out');
            xlim(xlimit);
            
            % Save the plot
            if i_ch == num_channels
                i_ch = 32;
            end
            if mod(i_ch-1,p)+1 == p
                xlabel(sprintf('%s       %s       %sHz       %s       Ch.%d-%d',group,region,frequency,mouse,i_ch-15,i_ch),'Interpreter','none');
                saveas(Fig1,fullfile(working_dir,'Figures_Indiv','PlotSpikeDelays_bySite',sprintf('%s_%s_%sHz_E%d.png',mouse,region,frequency,i_ch/16)),'png');
%                 pause;
                clf;
            end
        end        
    end
end
