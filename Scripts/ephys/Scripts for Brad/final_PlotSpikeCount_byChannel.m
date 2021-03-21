% Plot the number of the different neuronal responses by channel

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
namelist1 = {'INC','DEC','NC'};
namelist2 = {'low','high'};
linelist = {'-r','-b','-k'};

fileTag = '_filt';

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
All_Data = cell(length(Groups),length(Regions),2,16,3);
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_spikeTest.mat',region,frequency));
    if exist(Filename,'file')
        fprintf('%d - Loading: %s\n',i,Filename);
        load(Filename);
        
        [~,i_group] = ismember(group,Groups);
        [~,i_region] = ismember(region,Regions);
        if frequency == '20'
            i_frequency = 2;
        else
            i_frequency = 1;
        end
                
        % Get the spike counts
        fprintf('  Calculating spike counts...\n');
        Inc = sum(INC,3);
        Dec = sum(DEC,3);
        Nc = sum(NC,3);
        All_Data{i_group,1,i_frequency,1} = [All_Data{i_group,1,i_frequency,1};Inc(1,1:16)];
        All_Data{i_group,1,i_frequency,2} = [All_Data{i_group,1,i_frequency,2};Dec(1,1:16)];
        All_Data{i_group,1,i_frequency,3} = [All_Data{i_group,1,i_frequency,3};Nc(1,1:16)];
        All_Data{i_group,i_region,i_frequency,1} = [All_Data{i_group,i_region,i_frequency,1};Inc(1,17:32)];
        All_Data{i_group,i_region,i_frequency,2} = [All_Data{i_group,i_region,i_frequency,2};Dec(1,17:32)];
        All_Data{i_group,i_region,i_frequency,3} = [All_Data{i_group,i_region,i_frequency,3};Nc(1,17:32)];
    end
end

for i_freq = 1:2
    for i_response = 1:3
        clf;
        p = 1;
        for i_group = 1:size(All_Data,1)
            for i_region = 1:size(All_Data,2)
                subplot(5,5,p);
                clear Data;
                Data = All_Data{i_group,i_region,i_freq,i_response};
                if ~isempty(Data)
                    % Plot the results
                    plot([0 17],[0 0],'k','LineWidth',1);
                    hold on;
                    for i_mouse = 1:size(Data,1)
%                         plot([1:16],Data(i_mouse,:));
                        s = scatter([1:16],Data(i_mouse,:),'x');
%                         s.MarkerEdgeAlpha = 0.3;
                    end
                    
                    groupMean = mean(Data,1);
                    plot([1:16],groupMean,linelist{i_response},'LineWidth',1.5);
%                     scatter([1:16],groupMedian,'filled','or');
                    
                    set(gca,'XTick',[1:16]);
                    set(gca,'XTickLabel',[17:32],'fontsize',8);
                    set(gca,'YTick',[0:100]);
                    grid on;
                    xlim([0 17]);
                    
                    ylimit = get(gca,'YLim');
                    ylim([-0.5 ylimit(2)+0.5]);
                    
                    
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
        fprintf('Saving: %s\n',fullfile(working_dir,'Figures_Group',sprintf('PlotSpikeCount_%s_%s.png',namelist1{i_response},namelist2{i_freq})));
        saveas(Fig1,fullfile(working_dir,'Figures_Group',sprintf('PlotSpikeCount_%s_%s.png',namelist1{i_response},namelist2{i_freq})),'png');
%         pause;
    end
end