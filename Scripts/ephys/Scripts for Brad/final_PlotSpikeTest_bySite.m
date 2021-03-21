% Plot the firing rates per period per neuron for all neurons colored by
% response type by electrode site

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_spikeTest.mat',region,frequency));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        
        num_trials = size(Rates,1);
        num_channels = size(Rates,2);
        num_units = size(Rates,3);
        
        % Loop through and plot all of the neurons
        for i_ch = 1:num_channels
            
            p = 16;
            subplot(4,4,mod(i_ch-1,p)+1);
            hold on;
            
            for i_unit = 1:num_units
                if ~isempty(avg_ratesPerPeriod{1,i_ch,i_unit})
                    if INC(1,i_ch,i_unit) == 1
                        marker = '-r';
                    elseif DEC(1,i_ch,i_unit) == 1
                        marker = '-b';
                    elseif NC(1,i_ch,i_unit) == 1
                        marker = '--k';
                    end
                    errorbar([0.5 1.5 2.5],avg_ratesPerPeriod{1,i_ch,i_unit},sem_ratesPerPeriod{1,i_ch,i_unit},marker);
                    text(0.02,0.90,['Ch.' char(num2str(i_ch))],'Units','Normalized');
                end
            end
            
            % Format the plots
            set(gca,'XTick',[0.5 1.5 2.5]);
            set(gca,'XTickLabel',{'Pre','Stim','Post'});
            xlim([0 3]);
            ylabel('Hz');
            
            % Save the plot
            if i_ch == num_channels
                i_ch = 32;
            end
            if mod(i_ch-1,p)+1 == p
                xlabel(sprintf('%s       %s       %sHz       %s       Ch.%d-%d',group,region,frequency,mouse,i_ch-15,i_ch),'Interpreter','none');
                saveas(Fig1,fullfile(working_dir,'Figures_Indiv','PlotSpikeTest_bySite',sprintf('%s_%s_%sHz_E%d.png',mouse,region,frequency,i_ch/16)),'png');
%                 pause;
                clf;
            end
        end
    end
end