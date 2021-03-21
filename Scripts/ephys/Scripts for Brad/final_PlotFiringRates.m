% Plot the firing rates over time for all neurons grouped by type of change during stim

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

%%
color_area = {[173 0 0]./255,[0 0 173]./255,[0 173 0]./255,[0 0 0]./255};
color_line = {[215 0 0]./255,[0 0 215]./255,[0 215 0]./255,[128 128 128]./255};
names = {'Inc','Dec','NC','All'};
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
        
        % Loop through and separate out the firing rates for each electrode
        Elec_ch = {[1:16],[17:size(Rates,2)]};
        for i_elec = 1:2
            IncRates = [];
            DecRates = [];
            NCRates = [];
            AllRates = [];
            for i_channel = Elec_ch{i_elec}
                for i_unit = 1:size(Rates,3)
                    if INC(1,i_channel,i_unit) == 1
                        IncRates = [IncRates;avg_ratesPerBin{1,i_channel,i_unit}];
                        AllRates = [AllRates;avg_ratesPerBin{1,i_channel,i_unit}];
                    elseif DEC(1,i_channel,i_unit) == 1
                        DecRates = [DecRates;avg_ratesPerBin{1,i_channel,i_unit}];
                        AllRates = [AllRates;avg_ratesPerBin{1,i_channel,i_unit}];
                    elseif NC(1,i_channel,i_unit) == 1
                        NCRates = [NCRates;avg_ratesPerBin{1,i_channel,i_unit}];
                        AllRates = [AllRates;avg_ratesPerBin{1,i_channel,i_unit}];
                    end
                end
            end
            Rates_cell = {IncRates,DecRates,NCRates,AllRates};
            % Plot the firing rates
            OP.handle = figure(1);
            OP.alpha = 0.2;
            OP.line_width = 2;
            OP.error = 'sem';
            for i_rates = 1:length(Rates_cell)
                OP.color_area = color_area{i_rates};
                OP.color_line = color_line{i_rates};
                plot_areaerrorbar(Rates_cell{i_rates},OP);
                hold on;
                ylimits = ylim;
                ymax = ylimits(2);
                ylimits = [0 ymax*1.5];
                line([20 20],ylimits,'Color','k');
                line([40 40],ylimits,'Color','k');
                
                ylim(ylimits);
                xlim([0 60]);
                title({sprintf('%s:   %s   %s   %sHz',group,mouse,region,frequency),sprintf('Electrode %d - %s Units',i_elec,names{i_rates})},'Interpreter','none');
                ylabel('Firing Rate (Hz)');
                xlabel('seconds');
                
                saveas(Fig1,fullfile(working_dir,'Figures_Indiv','PlotFiringRates',sprintf('%s_%s_%sHz_E%d_%s.png',mouse,region,frequency,i_elec,names{i_rates})),'png');
%                 pause;
                clf;
            end
            clear Rates_cell;
        end
    end
end