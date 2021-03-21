
% Plot the stacked bar graph for the different neuronal responses

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Grad School\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

tag = '_5d';

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
hold on;
Data = cell(length(Groups),length(Regions),2);
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_spikeTest%s.mat',region,frequency,tag));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        
        % Calculate values for the stacked bar graphs
        E1_nINC = sum(sum(INC(:,1:16,:)));
        E2_nINC = sum(sum(INC(:,17:end,:)));
        E1_nDEC = sum(sum(DEC(:,1:16,:)));
        E2_nDEC = sum(sum(DEC(:,17:end,:)));
        E1_nNC = sum(sum(NC(:,1:16,:)));
        E2_nNC = sum(sum(NC(:,17:end,:)));
        E1_total = E1_nINC + E1_nDEC + E1_nNC;
        E2_total = E2_nINC + E2_nDEC + E2_nNC;
        
        [~,i_group] = ismember(group,Groups);
        [~,i_region] = ismember(region,Regions);
        if frequency == '20'
            i_frequency = 2;
        else
            i_frequency = 1;
        end
        Data{i_group,1,i_frequency} = [Data{i_group,1,i_frequency};E1_nINC,E1_nDEC,E1_nNC];
        Data{i_group,i_region,i_frequency} = [Data{i_group,i_region,i_frequency};E2_nINC,E2_nDEC,E2_nNC];
    end
end

p = 1;
for i_group = 1:size(Data,1)
    group = Groups{i_group};
    group = strsplit(group,'_');
    if group{2}(1) == 'C'
        tmp_freq = '6Hz';
    else
        tmp_freq = '10Hz';
    end
    for i_region = 1:size(Data,2)
        subplot(5,5,p);
        if ~isempty(Data{i_group,i_region,1})
            for i_frequency = 1:size(Data,3)
                totalInc(i_frequency) = sum(Data{i_group,i_region,i_frequency}(:,1));
                totalDec(i_frequency) = sum(Data{i_group,i_region,i_frequency}(:,2));
                totalNC(i_frequency) = sum(Data{i_group,i_region,i_frequency}(:,3));
            end
            totalAll = totalInc + totalDec + totalNC;
            pctInc = 100*totalInc./totalAll;
            pctDec = 100*totalDec./totalAll;
            pctNC = 100*totalNC./totalAll;
                        
            % Create the stacked bar graphs
            B = bar([pctInc(1),pctDec(1),pctNC(1);pctInc(2),pctDec(2),pctNC(2)],'stacked');
            B(1).FaceColor = 'r';
            B(2).FaceColor = 'b';
            B(3).FaceColor = [128 128 128]/255;
            set(gca,'XTick',[1 2]);
            set(gca,'TickDir','out'); box off;
            set(gca,'XTickLabel',{sprintf('%s (n=%d)',tmp_freq,totalAll(i_frequency)),'20Hz'});
            xlim([0.5 2.5]);
            set(gca,'YTick',[0:20:100]);
            ylim([0 100]);
            set(gca,'color','none');
            if i_region == 1
                ylabel(Groups{i_group},'Interpreter','none');
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
saveas(Fig1,fullfile(working_dir,'Figures_Group',sprintf('StackedBar_All%s.png',tag)),'png');