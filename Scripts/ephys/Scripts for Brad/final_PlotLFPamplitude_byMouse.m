% Plot the LFP amplitude for each mouse across all groups

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

Groups = {'Old ChAT','Old ChAT/APP','Young ChAT','Young VGAT','Young VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
namelist = {'low','high'};
comparelist = {'ChATs','CellTypes','All'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');

fprintf('Loading: LFPamplitude.mat\n');
Filename = fullfile(working_dir,'LFPamplitude.mat');

load(Filename);

%%

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
% %                     for i_mouse = 1:size(All_Data,5)
% %                         temp_data = squeeze(squeeze(All_Data(i_group,i_region,i_freq,:,i_mouse)));
% %                         if ~all(cellfun('isempty',temp_data))
% %                             temp_data = reshape(cell2mat(temp_data'),[size(temp_data{1},1),6,16-sum(cellfun('isempty',temp_data))]);
% %                             temp_data = temp_data(:,2,:)./temp_data(:,1,:);
% %                             temp_data = mean(squeeze(temp_data),2,'omitnan');
% %                             temp_data = temp_data - 1;
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
% %                     scatter(xlist,avglist,'k','filled');
% %                     errorbar(xlist,avglist,semlist,'k','LineStyle','none');
% %                     plot([i_group-0.2 i_group],[median(avglist) median(avglist)],'Color',[0.4940 0.1840 0.5560],'Linewidth',1.5);
% %                     plot([i_group-0.2 i_group+0.2],[mean(avglist) mean(avglist)],'Color',[0.4660 0.6740 0.1880],'Linewidth',1.5);
% %                     plot([i_group i_group+0.2],[median(avglist) median(avglist)],'Color',[0.4940 0.1840 0.5560],'Linewidth',1.5);
% %                 end
% %                 plot([0 6],[0 0],'-k');
% %                 ylimits = ylim;
% % %                 text(grouprange-.2,repmat(ylimits(2)*.05,1,3),'Median','HorizontalAlignment','center');
% % %                 text(grouprange+.2,repmat(ylimits(2)*.05,1,3),'Mean','HorizontalAlignment','center')
% %                 set(gca,'XTick',1:5);
% %                 set(gca,'TickDir','out');
% %                 set(gca,'XTickLabels',Groups);
% %                 xlim([0.5+d 3.5+d]);
% %                 title(sprintf('%s %s %s',comparelist{i_compare},Regions{i_region},namelist{i_freq}));
% %                 SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('PlotLFPamplitude_byMouse_%s_%s_%s.png',comparelist{i_compare},Regions{i_region},namelist{i_freq}));
% %                 fprintf('Saving: %s\n',SaveFilename);
% %                 saveas(Fig1,SaveFilename,'png');
% % %                 pause;
% %             end
% %         end
% %     end
% % end


dlist = [0,2,0];
labellist = {'6Hz','6Hz','6Hz','10Hz','10Hz'};
for i_compare = 3
    d = dlist(i_compare);
    grouprange = (1+d):(5+d);
    for i_region = 1:size(All_Data,2)
        clf;
        tt = 1;
        for i_group = grouprange
            if i_region ~= 3 || ~ismember(i_group,[1,2])
                avglist = [];
                semlist = [];
                for i_mouse = 1:size(All_Data,5)
                    temp_data = squeeze(squeeze(All_Data(i_group,i_region,1,:,i_mouse)));
                    if ~all(cellfun('isempty',temp_data))
                        temp_data = reshape(cell2mat(temp_data'),[size(temp_data{1},1),6,16-sum(cellfun('isempty',temp_data))]);
                        temp_data = temp_data(:,2,:)./temp_data(:,1,:);
                        temp_data = mean(squeeze(temp_data),2,'omitnan');
                        temp_data = temp_data - 1;
                        avg(1) = mean(temp_data);
                        sd = std(temp_data);
                        sem(1) = sd/sqrt(length(temp_data));
                        temp_data = squeeze(squeeze(All_Data(i_group,i_region,2,:,i_mouse)));
                        temp_data = reshape(cell2mat(temp_data'),[size(temp_data{1},1),6,16-sum(cellfun('isempty',temp_data))]);
                        temp_data = temp_data(:,2,:)./temp_data(:,1,:);
                        temp_data = mean(squeeze(temp_data),2,'omitnan');
                        temp_data = temp_data - 1;
                        avg(2) = mean(temp_data);
                        sd = std(temp_data);
                        sem(2) = sd/sqrt(length(temp_data));
                        p = errorbar([i_group-0.2 i_group+0.2],avg,sem,'--k');
                        hold on;
                        avglist(:,end+1) = avg;
                    end
                end
                plot([i_group-0.2 i_group+0.2],[median(avglist,2)],'Color',[0.8500 0.3250 0.0980],'Linewidth',2);
                plot([i_group-0.2 i_group+0.2],[mean(avglist,2)],'Color',[0.4660 0.6740 0.1880],'Linewidth',2);
                t = median(avglist,2);
                plot([i_group i_group+0.2],[mean(t) t(2)],'Color',[0.8500 0.3250 0.0980],'Linewidth',1.5);
                text((1+2*(tt-1))/10-1/25,0.08,labellist{i_group},'Units','normalized','HorizontalAlignment','center');
                text((1+2*(tt-1))/10+1/25,0.08,'20Hz','Units','normalized','HorizontalAlignment','center');
                text(.02,.07,'Mean','Color',[0.4660 0.6740 0.1880],'FontSize',12,'Units','normalized');
                text(.02,.04,'Median','Color',[0.8500 0.3250 0.0980],'FontSize',12,'Units','normalized');
            end
            tt = tt+1;
        end
        plot([0 6],[0 0],'-k');
        ylimits = ylim;
        set(gca,'XTick',1:5);
        set(gca,'TickDir','out');
        set(gca,'XTickLabels',Groups);
        ylabel('LFP Change (x100%)');
        xlim([0.5+d 5.5+d]);
        title(sprintf('%s LFP Change',Regions{i_region}));
        
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        left = outerpos(1) + ti(1)*2;
        bottom = outerpos(2) + ti(2)*2;
        ax_width = outerpos(3) - ti(1)*4 - ti(3);
        ax_height = outerpos(4) - ti(2)*2 - ti(4)*2;
        ax.Position = [left bottom ax_width ax_height];
        
        SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('PlotLFPamplitude_byMouse_%s_%s.png',Regions{i_region},comparelist{i_compare}));
        fprintf('Saving: %s\n',SaveFilename);
        saveas(Fig1,SaveFilename,'png');
%         pause;
    end
end