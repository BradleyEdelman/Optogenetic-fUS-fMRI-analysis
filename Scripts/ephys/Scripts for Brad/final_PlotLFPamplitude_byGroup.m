% Plot the LFP amplitude for each mouse at stimulation frequency and save them

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
namelist = {'low','high'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');

fprintf('Loading: LFPamplitude.mat\n');
Filename = fullfile(working_dir,'LFPamplitude.mat');

load(Filename);

for i_freq = 1:2
    p = 1;
    for i_group = 1:size(All_Data,1)
        for i_region = 1:size(All_Data,2)
            subplot(5,5,p);
            clear Data temp_data Data_avg Data_sem;
            Data = squeeze(All_Data(i_group,i_region,i_freq,:,:));
            if ~all(all(cellfun('isempty',Data)))
                
                % Rearrange the data and run calculations
                for i_mouse = 1:size(Data,2)
                    for i_channel = 1:size(Data,1)
                        temp_data = Data{i_channel,i_mouse};
                        if ~isempty(temp_data)
                            
                            amp_change = temp_data(:,2)./temp_data(:,1) - 1;
%                             amp_change = temp_data(:,1)*1000;
                            
                            avg = mean(amp_change);
%                             avg = median(amp_change);
                            sd = std(amp_change);
                            sem = sd / sqrt(length(amp_change));
                            Data_avg(i_mouse,i_channel) = avg;
                            Data_sem(i_mouse,i_channel) = sem;
                        end
                    end 
                end
                
                Data_avg(Data_avg == 0) = NaN;
                
                groupMedian = median(Data_avg,1,'omitnan');
                
                % Plot the results
                f = plot([0 17],[0 0],'k','LineWidth',2);
%                 f.Color = [0.4660 0.6740 0.1880];
%                 f.Color(4) = 0.5;

                hold on;
                for i_mouse = 1:size(Data_avg,1)
%                     plot([1:16],Data_avg(i_mouse,:));
%                     errorbar([1:16],Data_avg(i_mouse,:),Data_sem(i_mouse,:));
                    s = scatter([1:16],Data_avg(i_mouse,:),'ob');
                    s.MarkerEdgeAlpha = 0.3;
                end
                
                g = plot([1:16],groupMedian,'r','LineWidth',1.5);
                
                set(gca,'XTick',[1:16]);
                set(gca,'XTickLabel',[17:32],'fontsize',8);
                set(gca','TickDir','out');
                xlim([0 17]);
                
                ylimit = get(gca,'YLim');
                ylim([min(ylimit(1),-0.5) max(ylimit(2),1)]);
                
                
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
    SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('PlotLFPamplitude_groupMedian_%s.png',namelist{i_freq}));
    fprintf('Saving: %s\n',SaveFilename);
    saveas(Fig1,SaveFilename,'png');
%     pause;
    clf;
end