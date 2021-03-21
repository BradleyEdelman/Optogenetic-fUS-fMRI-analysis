% Plot the LFP traces averaged across trials by electrode site

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

fileTag = '_filt';

tsrange = 10000:50000; % milliseconds
xrange = [-inf inf];
% opt_yrange = [-0.7 0.7];

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        
        for xx = 1:size(LFP_mV_adj,1)
            for yy = 1:size(LFP_mV_adj,2)
                if isempty(LFP_mV_adj{xx,yy})
                    LFP_mV_adj{xx,yy} = zeros(1,60000);
                end
            end
        end
        
        TS_cell = LFP_mV_adj;
        
        % Cut off the last few points of the timeseries due to some being 1ms longer or shorter
        TS_cell = cellfun(@(x) x(1:59999),TS_cell,'UniformOutput',false);
        TS_length = size(TS_cell{1,1},1);
        num_trials = size(TS_cell,1);
        num_channels = size(TS_cell,2);
        
        % Plot/Save the averages
        fprintf('Plotting/Saving: LFP Averages across Trials...\n');
        for i_ch = 1:num_channels
            TS = reshape(cell2mat(TS_cell(:,i_ch)),TS_length,[])';
            TS_avg = mean(TS,1);
            Laser_TS = Laser_TS_adj{1,1};
            LFP_TS = LFP_TS_adj{1,1};
            
            % Plot the LFP traces
            p = 16;
            subplot(8,2,mod(i_ch-1,p)+1);
            hold on;
            
            plot(LFP_TS(tsrange),TS_avg(tsrange),'k');
            if exist('opt_yrange','var')
                yrange = opt_yrange;
            else
                yrange = ylim;
            end
            scatter(Laser_TS,(yrange(1)+0.1)*ones(length(Laser_TS),1),'.b');
            text(LFP_TS(20000)-2,yrange(2)-0.15,['Ch.' char(num2str(i_ch))]);
            
            % Format the plots
            pos = get(gca, 'Position');
            pos(1) = 0.02+0.500*(mod(i_ch-1,2));
            pos(3) = 0.470;
            pos(4) = 0.08;
            set(gca, 'Position', pos);
            set(gca,'XTick',[0:10:60]);
            if exist('opt_yrange','var')
                set(gca,'YTick',[-1:0.5:1]);
            end
            xlim(xrange);
            ylim(yrange);
            
            % Save the plot
            if mod(i_ch-1,p)+1 == p
                xlabel(sprintf('%s       %s       %sHz       %s       Ch.%d-%d',group,region,frequency,mouse,i_ch-15,i_ch),'Interpreter','none');
                saveas(Fig1,fullfile(working_dir,'Figures_Indiv','PlotLFPaverages_bySite',sprintf('%s_%s_%sHz_E%d.png',mouse,region,frequency,i_ch/16)),'png');
%                 pause;
                clf;
            end
        end
    end
end