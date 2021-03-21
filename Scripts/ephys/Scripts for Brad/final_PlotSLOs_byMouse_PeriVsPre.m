% Plot the SLO count during vs. before stim for all groups

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Groups2 = {'Old ChAT','Old ChAT/APP','Young ChAT','Young VGAT','Young VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
Frequencies = {'low','high'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
All_Data = cell(length(Groups),length(Regions),2,2);
mouseTracker = {};
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_countSLOs_PrePeriPost.mat',region,frequency));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        
        [~,i_group] = ismember(group,Groups);
        [~,i_region] = ismember(region,Regions);
        if frequency == '20'
            i_frequency = 2;
        else
            i_frequency = 1;
        end
        
        i_mouse = find(cellfun(@(x) isequal([mouse region],x),mouseTracker),1);
        if i == 1
            i_mouse = 1;
        end
        if isempty(i_mouse)
            mouseTracker{end+1} = [mouse region];
            i_mouse = length(mouseTracker);
        end
        All_Data{i_group,1,i_frequency,i_mouse} = {SLOcount_pre(:,1:16),SLOcount_peri(:,1:16),SLOcount_post(:,1:16)};
        All_Data{i_group,i_region,i_frequency,i_mouse} = {SLOcount_pre(:,17:32),SLOcount_peri(:,17:32),SLOcount_post(:,17:32)};
    end
end

%%
for i_region = 1:size(All_Data,2)
    for i_freq = 1:size(All_Data,3)
        clf;
        tt = 1;
        for i_group = 1:size(All_Data,1)
            avglist = [];
            semlist = [];
            for i_mouse = 1:size(All_Data,4)
                if ~isempty(All_Data{i_group,i_region,i_freq,i_mouse})
                    
                    temp_data = All_Data{i_group,i_region,i_freq,i_mouse}{1};
                    temp_data = mean(temp_data,2,'omitnan');
                    avg(1) = mean(temp_data);
                    sd = std(temp_data);
                    sem(1) = sd/sqrt(length(temp_data));
                    
                    temp_data = All_Data{i_group,i_region,i_freq,i_mouse}{2};
                    temp_data = mean(temp_data,2,'omitnan');
                    avg(2) = mean(temp_data);
                    sd = std(temp_data);
                    sem(2) = sd/sqrt(length(temp_data));
                   
                    p = errorbar([i_group-0.2 i_group+0.2],avg,sem,'--k');
                    hold on;
                    avglist(:,end+1) = avg;
                end
            end
            if ~isempty(avglist)
                plot([i_group-0.2 i_group+0.2],[median(avglist,2)],'Color',[0.8500 0.3250 0.0980],'Linewidth',2);
                plot([i_group-0.2 i_group+0.2],[mean(avglist,2)],'Color',[0.4660 0.6740 0.1880],'Linewidth',2);
                t = median(avglist,2);
                plot([i_group i_group+0.2],[mean(t) t(2)],'Color',[0.8500 0.3250 0.0980],'Linewidth',1.5);
                text((1+2*(tt-1))/10-1/25,0.08,'Baseline','Units','normalized','HorizontalAlignment','center');
                text((1+2*(tt-1))/10+1/25,0.08,'Stimulation','Units','normalized','HorizontalAlignment','center');
            end
            tt = tt+1;
        end
        text(.01,.06,'Mean','Color',[0.4660 0.6740 0.1880],'FontSize',12,'Units','normalized');
        text(.01,.03,'Median','Color',[0.8500 0.3250 0.0980],'FontSize',12,'Units','normalized');
        plot([0 6],[0 0],'-k');
        set(gca,'XTick',1:5);
        set(gca,'YTick',-0.2:0.2:2.0);
        set(gca,'TickDir','out');
        set(gca,'XTickLabels',Groups2);
        ylabel('# of SLOs');
        xlim([0.5 5.5]);
        ylim([-0.2 2.0]);
        title(sprintf('%s SLOs %s freq',Regions{i_region},Frequencies{i_freq}));

        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        left = outerpos(1) + ti(1)*2;
        bottom = outerpos(2) + ti(2)*2;
        ax_width = outerpos(3) - ti(1)*4 - ti(3);
        ax_height = outerpos(4) - ti(2)*2 - ti(4)*2;
        ax.Position = [left bottom ax_width ax_height];

        SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('PlotSLOs_byMouse_PeriVsPre_%s_%s.png',Regions{i_region},Frequencies{i_freq}));
        fprintf('Saving: %s\n',SaveFilename);
        saveas(Fig1,SaveFilename,'png');
%         pause;
    end
end