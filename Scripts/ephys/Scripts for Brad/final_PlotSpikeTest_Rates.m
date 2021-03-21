% Plot the distribution of firing rates, separated by response type and
% stimulation period

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Grad School\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
Frequencies = {'low','high','both'};
Responses = {'INC','DEC','NC','ALL'};
Periods = {'Pre','Peri','Post'};

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','on');
hold on;
All_Data = cell(length(Groups),length(Regions),3,3,3,1);
n_unit = 1;
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
                        i_response = 1;
                    elseif DEC(1,i_ch,i_unit) == 1
                        i_response = 2;
                    elseif NC(1,i_ch,i_unit) == 1
                        i_response = 3;
                    end
                    for temp_freq = [i_frequency,3]
                        for temp_response = [i_response,4]
                            for temp_period = 1:3
                                if i_ch <= 16
                                    All_Data{i_group,1,temp_freq,temp_response,temp_period,n_unit} = avg_ratesPerPeriod{1,i_ch,i_unit}(temp_period);
                                else
                                    All_Data{i_group,i_region,temp_freq,temp_response,temp_period,n_unit} = avg_ratesPerPeriod{1,i_ch,i_unit}(temp_period);
                                end
                            end
                        end
                    end
                    n_unit = n_unit + 1;
                end
            end
        end
    end
end

for i_period = 1:size(All_Data,5)
    for i_response = 1:size(All_Data,4)
        for i_freq = 1:size(All_Data,3)
            p = 1;
            for i_group = 1:size(All_Data,1)
                for i_region = 1:size(All_Data,2)
                    temp_data = [];
                    for i_unit = 1:size(All_Data,6)                        
                        if ~isempty(All_Data{i_group,i_region,i_freq,i_response,i_period,i_unit})
                            temp_data = [temp_data , All_Data{i_group,i_region,i_freq,i_response,i_period,i_unit}];
                        end
                    end
                    figs(i_group,i_region) = subplot(5,5,p);
                    if ~isempty(temp_data)
                        % Plot the data
                        histogram(temp_data,'BinWidth',1);
                        hold on;
                        xlabel(sprintf('Mean:%.1f  Med:%.1f  n=%d',mean(temp_data),median(temp_data),length(temp_data)));


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
            for i_region = 1:size(All_Data,2)
                linkaxes(figs(:,i_region),'x');
            end
            % Save the plot
            SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('SpikeTest_Rates_%s_%s_%s.png',Periods{i_period},Responses{i_response},Frequencies{i_freq}));
            fprintf('Saving: %s\n',SaveFilename);
%             saveas(Fig1,SaveFilename,'png');
            pause;
            clf;
        end
    end
end