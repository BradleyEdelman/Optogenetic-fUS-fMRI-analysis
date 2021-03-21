% Plot the LFP peaks

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

fprintf('Loading: LFPpeaks.mat\n');
Filename = fullfile(working_dir,'LFPpeaks.mat');
load(Filename);

for i_ch = 1:size(All_Data,4)
    for i_frequency = 1:size(All_Data,3)
        p = 1;
        for i_group = 1:size(All_Data,1)
            for i_region = 2:size(All_Data,2)
                subplot(5,4,p);
                if ~isempty(All_Data{i_group,i_region,i_frequency,i_ch})
                    
                    pts = All_Data{i_group,i_region,i_frequency,i_ch};
                    ylimits = [0 0.001];
                    for kk = 1:size(pts,1)
                        if sum(pts(kk,:)) ~= 0
                            g = plot(pts(kk,401:600));
                            hold on;
                            g.Color(4) = 0.7;
                            new_ylimits = ylim;
                            ylimits(1) = min(ylimits(1),new_ylimits(1));
                            ylimits(2) = max(ylimits(2),new_ylimits(2));
                        end
                    end
                    verts = [100,105,110,115,120,125,130,135,140,145];
                    for i_v = 1:length(verts)
                        p1 = line([verts(i_v) verts(i_v)],[ylimits(1) ylimits(2)],'Color','k');
                        p1.Color(4) = 0.3;
                    end
                    ylim(ylimits);
                    xlim([95 150]);
                    
                else
                    set(gca,'Visible','off');
                    ax = gca;
                    ax.Title.Visible = 'on';
                end
                if i_region == 2
                        ylabel(Groups{i_group},'Interpreter','none');
                        ax = gca;
                        ax.YLabel.Visible = 'on';
                end
                set(gca,'XTick',[0:5:200]);
                set(gca,'XTickLabel',[-100:5:100]);
                if p < 5
                    if i_ch <= 16
                        title(sprintf('%s Ch.%d (BF)',Regions{i_region},i_ch));
                    else
                        title(sprintf('%s Ch.%d',Regions{i_region},i_ch));
                    end
                end
                p = p+1;
            end
        end
        % Save the plot
        fprintf('Saving: %s',fullfile(working_dir,'Figures_Group','PlotLFPpeaks_Zoomed',sprintf('LFPpeaks_Ch%d_%s_zoomed.png\n',i_ch,namelist{i_frequency})));
        saveas(Fig1,fullfile(working_dir,'Figures_Group','PlotLFPpeaks_Zoomed',sprintf('LFPpeaks_Ch%d_%s_zoomed.png',i_ch,namelist{i_frequency})),'png');
%         pause;
        clf;
    end
end