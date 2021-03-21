% Plot the SLO count for each mouse across all groups

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Groups2 = {'Old ChAT','Old ChAT/APP','Young ChAT','Young VGAT','Young VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
namelist = {'low','high'};
comparelist = {'ChATs','CellTypes'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
All_Data = cell(length(Groups),length(Regions),2,2);
mouseTracker = {};
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_countSLOs.mat',region,frequency));
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
        All_Data{i_group,1,i_frequency,i_mouse} = SLOcount(:,1:16);
        All_Data{i_group,i_region,i_frequency,i_mouse} = SLOcount(:,17:32);
    end
end

% % dlist = [0,2];
% % for i_compare = 1:2
% %     d = dlist(i_compare);
% %     grouprange = (1+d):(3+d);
% %     for i_region = 1:size(All_Data,2)
% %         if i_compare == 2 || i_region ~= 3
% %             for i_freq = 1:2
% %                 clf;
% %                 for i_group = grouprange
% %                     avglist = [];
% %                     semlist = [];
% %                     for i_mouse = 1:size(All_Data,4)
% %                         temp_data = All_Data{i_group,i_region,i_freq,i_mouse};
% %                         if ~isempty(temp_data)
% %                             temp_data = mean(temp_data,2,'omitnan');
% %                             avg = mean(temp_data);
% %                             sd = std(temp_data);
% %                             sem = sd/sqrt(length(temp_data));
% %                             avglist(end+1) = avg;
% %                             semlist(end+1) = sem;
% %                         end
% %                     end
% %                     [avglist,I] = sort(avglist);
% %                     semlist = semlist(I);
% %                     xlist = i_group+(rand(1,length(avglist))-0.5)/5;
% %                     hold on;
% % %                     bar(i_group-.2,median(avglist),0.4,'EdgeColor','r','FaceColor','none');
% % %                     bar(i_group+.2,mean(avglist),0.4,'EdgeColor','b','FaceColor','none');
% % %                     text(i_group-0.2,median(avglist),'Median','HorizontalAlignment','center')
% % %                     text(i_group+0.2,mean(avglist),'Mean','HorizontalAlignment','center')
% %                     scatter(xlist,avglist,'k','filled');
% %                     errorbar(xlist,avglist,semlist,'k','LineStyle','none');
% %                     plot([i_group-0.2 i_group],[median(avglist) median(avglist)],'Color',[0.4940 0.1840 0.5560],'Linewidth',1.5);
% %                     plot([i_group-0.2 i_group+0.2],[mean(avglist) mean(avglist)],'Color',[0.4660 0.6740 0.1880],'Linewidth',1.5);
% %                     plot([i_group i_group+0.2],[median(avglist) median(avglist)],'Color',[0.4940 0.1840 0.5560],'Linewidth',1.5);
% %                 end
% %                 plot([0 6],[0 0],'-k');
% % %                 text(grouprange-.2,repmat(ylimits(2)*.05,1,3),'Median','HorizontalAlignment','center');
% % %                 text(grouprange+.2,repmat(ylimits(2)*.05,1,3),'Mean','HorizontalAlignment','center')
% %                 set(gca,'XTick',1:5);
% %                 set(gca,'YTick',-0.2:0.2:2.0);
% %                 set(gca,'TickDir','out');
% %                 set(gca,'XTickLabels',Groups2);
% %                 xlim([0.5+d 3.5+d]);
% %                 ylim([-0.2 2.0]);
% %                 title(sprintf('%s %s %s',comparelist{i_compare},Regions{i_region},namelist{i_freq}));
% %                 SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('PlotSLOs_byMouse_%s_%s_%s.png',comparelist{i_compare},Regions{i_region},namelist{i_freq}));
% %                 fprintf('Saving: %s\n',SaveFilename);
% %                 saveas(Fig1,SaveFilename,'png');
% % %                 pause;
% %             end
% %         end
% %     end
% % end

dlist = [0,2];
labellist = {'6Hz','6Hz','6Hz','10Hz','10Hz'};
for i_region = 1:size(All_Data,2)
    clf;
    tt = 1;
    for i_group = 1:5
        avglist = [];
        semlist = [];
        for i_mouse = 1:size(All_Data,4)
            temp_data = All_Data{i_group,i_region,1,i_mouse};
            if ~isempty(temp_data)
                temp_data = mean(temp_data,2,'omitnan');
                avg(1) = mean(temp_data);
                sd = std(temp_data);
                sem(1) = sd/sqrt(length(temp_data));
                temp_data = All_Data{i_group,i_region,2,i_mouse};
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
            text((1+2*(tt-1))/10-1/25,0.08,labellist{i_group},'Units','normalized','HorizontalAlignment','center');
            text((1+2*(tt-1))/10+1/25,0.08,'20Hz','Units','normalized','HorizontalAlignment','center');
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
    title(sprintf('%s SLOs',Regions{i_region}));
    
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1)*2;
    bottom = outerpos(2) + ti(2)*2;
    ax_width = outerpos(3) - ti(1)*4 - ti(3);
    ax_height = outerpos(4) - ti(2)*2 - ti(4)*2;
    ax.Position = [left bottom ax_width ax_height];
    
    SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('PlotSLOs_byMouse_%s.png',Regions{i_region}));
    fprintf('Saving: %s\n',SaveFilename);
    saveas(Fig1,SaveFilename,'png');
    %             pause;
end
