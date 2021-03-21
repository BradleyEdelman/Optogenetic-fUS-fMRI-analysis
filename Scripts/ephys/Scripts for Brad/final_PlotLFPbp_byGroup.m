% Plot the LFP band power for each mouse at stimulation frequency and save them

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
namelist = {'low','high'};

fileTag = '_filt';

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');

fprintf('Loading: LFPbp.mat\n');
Filename = fullfile(working_dir,'LFPbp.mat');

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
                            
                            bp_change = temp_data(:,2)./temp_data(:,1);
%                             bp_change = temp_data(:,2)*1000000;
                            
                            avg = mean(bp_change);
%                             avg = median(bp_change);
                            sd = std(bp_change);
                            sem = sd / sqrt(length(bp_change));
                            Data_avg(i_mouse,i_channel) = avg;
                            Data_sem(i_mouse,i_channel) = sem;
                        end
                    end 
                end
                
                % Plot the results
                hold on;
                for i_mouse = 1:size(Data_avg,1)
%                     plot([1:16],Data_avg(i_mouse,:));
                    errorbar([1:16],Data_avg(i_mouse,:),Data_sem(i_mouse,:));
                end
                
                set(gca,'XTick',[1:16]);
                set(gca,'XTickLabel',[17:32],'fontsize',8);
                xlim([0 17]);
                
%                 set(gca,'YScale','log');
%                 ylim([0 1000]);

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
    fprintf('Saving: %s\n',fullfile(working_dir,'Figures_Group',sprintf('PlotLFPbp_%s.png',namelist{i_freq})));
    saveas(Fig1,fullfile(working_dir,'Figures_Group',sprintf('PlotLFPbp_%s.png',namelist{i_freq})),'png');
%     pause;
    clf;
end