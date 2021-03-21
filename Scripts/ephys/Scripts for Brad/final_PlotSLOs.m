% Plot the SLO counts from the different animals

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
Data_AVG = cell(length(Groups),length(Regions),2);
Data_SEM = Data_AVG;
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_countSLOs.mat',region,frequency));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        
        % Calculate mean and SEM for SLOs        
        num_channels = size(SLOcount,2);
        num_trials = size(SLOcount,1);
        
        avg_nSLO = mean(SLOcount,1);
        sem_nSLO = std(SLOcount) / sqrt(num_trials);
        
        [~,i_group] = ismember(group,Groups);
        [~,i_region] = ismember(region,Regions);
        if frequency == '20'
            i_frequency = 2;
        else
            i_frequency = 1;
        end
        
        Data_AVG{i_group,1,i_frequency} = [Data_AVG{i_group,1,i_frequency};avg_nSLO(:,1:16)];
        Data_AVG{i_group,i_region,i_frequency} = [Data_AVG{i_group,i_region,i_frequency};avg_nSLO(:,17:32)];
        
        Data_SEM{i_group,1,i_frequency} = [Data_SEM{i_group,1,i_frequency};sem_nSLO(:,1:16)];
        Data_SEM{i_group,i_region,i_frequency} = [Data_SEM{i_group,i_region,i_frequency};sem_nSLO(:,17:32)];
    end
end

namelist = {'low','high'};
for i_freq = 1:2
    clf;
    p = 1;
    for i_group = 1:size(Data_AVG,1)
        for i_region = 1:size(Data_AVG,2)
            subplot(5,5,p);
            if ~isempty(Data_AVG{i_group,i_region,i_freq})
                
                % Plot the data
                current_avg = Data_AVG{i_group,i_region,i_freq};
                current_sem = Data_SEM{i_group,i_region,i_freq};
                hold on;                
                plot([0 17],[0 0],'k','LineWidth',1);
                for i_mouse = 1:size(current_avg,1)
%                     plot([1:16],current_avg(i_mouse,:));
                    errorbar([1:16],current_avg(i_mouse,:),current_sem(i_mouse,:));
%                     s = scatter([1:16],current_avg(i_mouse,:),'ob');
%                     s.MarkerEdgeAlpha = 0.5;
                end
                
%                 groupMedian = median(current_avg,1,'omitnan');
%                 g = plot([1:16],groupMedian,'r','LineWidth',1.5);
                
                xlabel(sprintf('Avg: %.2f SLOs',mean(groupMedian)));
                set(gca,'XTick',[1:16]);
                set(gca,'XTickLabel',[17:32],'fontsize',8);
                xlim([0 17]);
                set(gca,'YTick',[0:0.5:5]);               
                ylim([-0.5 2.5]);
                
                if i_region == 1
                    ylabel(Groups{i_group},'Interpreter','none');
                    set(gca,'XTickLabel',[1:16],'fontsize',8);
                end
            else
                set(gca,'Visible','off');
                ax = gca;
                ax.Title.Visible = 'on';
            end
            if p <= 5
                title(Regions{i_region});
            end
            p = p+1;
        end
    end
    % Save the plot
    SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('PlotSLOs_%s.png',namelist{i_freq}));
    fprintf('Saving: %s\n',SaveFilename);
    saveas(Fig1,SaveFilename,'png');
%     pause;
end