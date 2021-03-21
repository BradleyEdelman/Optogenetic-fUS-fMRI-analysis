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
                        i_region_tmp = 1;
                    else
                        i_region_tmp = i_region;
                    end
%                     fprintf('i:%d g:%d r:%d f:%d c:%d\n',i,i_group,i_region_tmp,i_frequency,i_ch);
                    i_ch_mod = (mod(i_ch-1,16)+1);
                    All_Data{i_group,i_region_tmp,i_frequency,i_ch_mod,n_unit} = {avg_ratesPerPeriod{1,i_ch,i_unit},sem_ratesPerPeriod{1,i_ch,i_unit},marker};
                    All_Data{i_group,i_region_tmp,3,i_ch_mod,n_unit} = {avg_ratesPerPeriod{1,i_ch,i_unit},sem_ratesPerPeriod{1,i_ch,i_unit},marker};
                    n_unit = n_unit + 1;
                end
            end
        end
    end
end

for i_freq = 1:3
    fprintf('f: %d\n',i_freq);
    for i_group = 1:size(All_Data,1)
        fprintf(' g: %d\n',i_group);
        for i_region = 1:size(All_Data,2)
            fprintf('  r: %d\n',i_region);
            if i_region > 1
                k = 16;
            else
                k = 0;
            end
            save_flag = 0;
            if any(any(~cellfun('isempty',All_Data(i_group,i_region,i_freq,:,:))))
                save_flag = 1;
                p = 1;
                for i_ch = 1:16
                    subplot(4,4,p);
                    nINC = 0;
                    nDEC = 0;
                    nNC = 0;
                    fprintf('   c: %d\n',i_ch);
                    % Plot the data
                    for i_unit = 1:size(All_Data,5)
                        if ~isempty(All_Data{i_group,i_region,i_freq,i_ch,i_unit})                           
                            current_unit = All_Data{i_group,i_region,i_freq,i_ch,i_unit};
                            
                            if isequal(current_unit{3},'-r')
                                nINC = nINC + 1;
                            elseif isequal(current_unit{3},'-b')
                                nDEC = nDEC + 1;
                            else
                                nNC = nNC + 1;
                            end
                            
                            errorbar([0.5 1.5 2.5],current_unit{1},current_unit{2},current_unit{3});
                            hold on;
                        end
                    end
                    
                    % Format the plots
                    set(gca,'XTick',[0.5 1.5 2.5]);
                    set(gca,'XTickLabel',{'Pre','Stim','Post'});
                    title(sprintf('Ch.%d',i_ch+k));
                    xlim([0 3]);
                    
                    nTotal = nINC + nDEC + nNC;
                    text(0.87,0.8,sprintf('I: %d',nINC),'Units','normalized');
                    text(0.87,0.6,sprintf('D: %d',nDEC),'Units','normalized');
                    text(0.87,0.4,sprintf('N: %d',nNC),'Units','normalized');
                    text(0.87,0.2,sprintf('T: %d',nTotal),'Units','normalized');
                    
                    p = p+1;
                end
            end
            % Save the plot
            if save_flag == 1
                xlabel(sprintf('%s       %s       %s_freq',Groups{i_group},Regions{i_region},namelist{i_freq}),'Interpreter','none');
                fprintf('Saving: %s',fullfile(working_dir,'Figures_Group','PlotSpikeTest_byChannel',sprintf('SpikeTest_%s_%s_%s.png\n',Groups{i_group},Regions{i_region},namelist{i_freq})));
                saveas(Fig1,fullfile(working_dir,'Figures_Group','PlotSpikeTest_byChannel',sprintf('SpikeTest_%s_%s_%s.png',Groups{i_group},Regions{i_region},namelist{i_freq})),'png');
%                 pause;
            end
            clf;
        end
    end
end