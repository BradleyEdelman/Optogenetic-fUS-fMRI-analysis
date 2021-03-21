% Get the LFP PSD for each mouse and save them

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
namelist = {'low','high'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
climit = [-135 -60];
ylimit = [4 30];

for i_group = 1:length(Groups)
    group = Groups{i_group};
    fprintf('Loading: LFPpsd_%s.mat\n\n',group);
    load(fullfile(working_dir,sprintf('LFPpsd_%s.mat',group)));
    for i_region = 1:length(Regions)
        region = Regions{i_region};
        if i_region == 1
            d = 0;
        else
            d = 16;
        end
        for i_frequency = 1:size(Data,3)
            if ~isempty(Data{1,i_region,1,1})
                
                fprintf('Plotting/Saving: %s %s %s\n',group,region,namelist{i_frequency});
                for i_channel = 1:16
                    subplot(4,4,i_channel);
                    
                    temp_data = mean(Data{1,i_region,i_frequency,i_channel},3);
                    imagesc(mag2db(temp_data)); axis xy
                    hold on;
                    text(0.02,0.90,['Ch.' char(num2str(i_channel+d))],'Units','Normalized','FontWeight','bold','Color','w');
                    % Format the plots
                    set(gca,'XTick',[0:99:297]);
                    set(gca,'XTickLabel',[0 20 40 60]);
                    ax=gca;
                    ax.CLim = climit;
                    ylim(ylimit);
                end

                % Save the plot
                xlabel(sprintf('%s       %s       %s',group,region,namelist{i_frequency}),'Interpreter','none');
                saveas(Fig1,fullfile(working_dir,'Figures_Group','PlotLFPpsd_byGroup',sprintf('LFPpsd_%s_%s_%s.png',region,group,namelist{i_frequency})),'png');
%                 pause;
                clf;
            end
        end
    end
end