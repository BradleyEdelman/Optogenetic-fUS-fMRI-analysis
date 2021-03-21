% Plot the firing rates per period per neuron for all neurons colored by
% response type, separate graphs for each geno and frequency

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
hold on;
All_Data = cell(length(Groups),length(Regions),3,1);
n_unit = 1;
namelist = {'low','high','both'};
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_spikeTest.mat',region,frequency));
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
        
        % Get the relevant data from the spikeTest files
        for i_ch = 1:size(avg_ratesPerPeriod,2)
            for i_unit = 1:size(avg_ratesPerPeriod,3)
                if ~isempty(avg_ratesPerPeriod{1,i_ch,i_unit})
                    if INC(1,i_ch,i_unit) == 1
                        marker = '-r';
                    elseif DEC(1,i_ch,i_unit) == 1
                        marker = '-b';
                    elseif NC(1,i_ch,i_unit) == 1
                        marker = '--k';
                    end
                    if i_ch <= 16
                        All_Data{i_group,1,i_frequency,n_unit} = {avg_ratesPerPeriod{1,i_ch,i_unit},sem_ratesPerPeriod{1,i_ch,i_unit},marker};
                        All_Data{i_group,1,3,n_unit} = {avg_ratesPerPeriod{1,i_ch,i_unit},sem_ratesPerPeriod{1,i_ch,i_unit},marker};
                    else
                        All_Data{i_group,i_region,i_frequency,n_unit} = {avg_ratesPerPeriod{1,i_ch,i_unit},sem_ratesPerPeriod{1,i_ch,i_unit},marker};
                        All_Data{i_group,i_region,3,n_unit} = {avg_ratesPerPeriod{1,i_ch,i_unit},sem_ratesPerPeriod{1,i_ch,i_unit},marker};
                    end
                    n_unit = n_unit + 1;
                end
            end
        end
    end
end

for i_freq = 1:3
    fprintf('f: %d\n',i_freq);
    for i_group = 1:size(All_Data,1)
        fprintf('g: %d\n',i_group);
        p = 1;
        for i_region = 1:size(All_Data,2)
            fprintf('r: %d\n',i_region);
            subplot(1,5,p);
            if any(~cellfun('isempty',All_Data(i_group,i_region,i_freq,:)))
                
                % Plot the data
                for i_unit = 1:size(All_Data,4)
                    if ~isempty(All_Data{i_group,i_region,i_freq,i_unit})
                        current_unit = All_Data{i_group,i_region,i_freq,i_unit};
                        errorbar([0.5 1.5 2.5],current_unit{1},current_unit{2},current_unit{3});
                        hold on;
                    end
                end
                % Format the plots
                pos = get(gca, 'Position');
                pos(2) = 0.05;
                pos(4) = 0.90;
                set(gca, 'Position', pos);
                set(gca,'XTick',[0.5 1.5 2.5]);
                set(gca,'XTickLabel',{'Pre','Stim','Post'});
                xlim([0 3]);
                
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
        % Save the plot
        fprintf('Saving: %s',fullfile(working_dir,'Figures_Group',sprintf('SpikeTest_%s_%s.png\n',Groups{i_group},namelist{i_freq})));
        saveas(Fig1,fullfile(working_dir,'Figures_Group',sprintf('SpikeTest_%s_%s.png',Groups{i_group},namelist{i_freq})),'png');
%         pause;
        clf;
    end
end

% for i_freq = 1:3
%     fprintf('f: %d\n',i_freq);
%     for i_region = 1:size(All_Data,2)
%         fprintf('r: %d\n',i_region);
%         p = 1;
%         for i_group = 1:size(All_Data,1)
%             fprintf('g: %d\n',i_group);
%             subplot(1,5,p);
%             if any(~cellfun('isempty',All_Data(i_group,i_region,i_freq,:)))
%                 
%                 % Plot the data
%                 for i_unit = 1:size(All_Data,4)
%                     if ~isempty(All_Data{i_group,i_region,i_freq,i_unit})
%                         current_unit = All_Data{i_group,i_region,i_freq,i_unit};
%                        
% %                         %Normalize
% %                         current_unit{1} = current_unit{1}/current_unit{1}(1);
% %                         current_unit{2} = current_unit{2}/current_unit{1}(1);
% %                         plot([0.5 1.5 2.5],current_unit{1},current_unit{3});
% 
%                         errorbar([0.5 1.5 2.5],current_unit{1},current_unit{2},current_unit{3});
%                         hold on;
%                     end
%                 end
%                 % Format the plots
%                 pos = get(gca, 'Position');
%                 pos(2) = 0.05;
%                 pos(4) = 0.90;
%                 set(gca, 'Position', pos);
%                 set(gca,'XTick',[0.5 1.5 2.5]);
%                 set(gca,'XTickLabel',{'Pre','Stim','Post'});
%                 xlim([0 3]);
%                 
%             else
%                 set(gca,'Visible','off');
%                 ax = gca;
%                 ax.Title.Visible = 'on';
%             end
%             if p <= 5
%                 title(Groups{i_group},'Interpreter','none');
%             end
%             p = p+1;
%         end
%         % Save the plot
%         fprintf('Saving: %s',fullfile(working_dir,'Figures_Group',sprintf('SpikeTest_%s_%s.png\n',Regions{i_region},namelist{i_freq})));
%         saveas(Fig1,fullfile(working_dir,'Figures_Group',sprintf('SpikeTest_%s_%s.png',Regions{i_region},namelist{i_freq})),'png');
% %         pause;
%         clf;
%     end
% end