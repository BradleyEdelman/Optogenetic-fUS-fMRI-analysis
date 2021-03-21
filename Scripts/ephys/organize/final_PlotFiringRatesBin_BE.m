function final_PlotFiringRatesBin_BE(storage,base_fold,slash)
% Plot the firing rates over time for all neurons grouped by type of change during stim

%%
color_area = {[173 0 0]./255,[0 0 173]./255,[0 173 0]./255,[0 0 0]./255};
color_line = {[215 0 0]./255,[0 0 215]./255,[0 215 0]./255,[128 128 128]./255};
names = {'Inc','Dec','NC','All'};
% Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');

stim = {'0_1' '0_5' '1_0'};
IncRates = cell(size(base_fold,1),size(stim,2),4);
DecRates = cell(size(base_fold,1),size(stim,2),4);
NCRates = cell(size(base_fold,1),size(stim,2),4);
AllRates = cell(size(base_fold,1),size(stim,2),4);
for i_mouse = 1:size(base_fold,1)
    
    for i_stim = 1:size(stim,2)
        
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '_spike_test.mat'];
        
        if exist(data_file,'file')
            load(data_file);
        
            % Loop through and separate out the firing rates for each region
            region_idx = {1:8;9:16;17:24;25:32};
            for i_region = 1:4
                for i_channel = 1:size(region_idx{i_region,:},2)
                    for i_unit = 1:size(Rates,3)
                        if INC(1,region_idx{i_region,:}(i_channel),i_unit) == 1
                            IncRates{i_mouse,i_stim,i_region} = [IncRates{i_mouse,i_stim,i_region};avg_ratesPerBin{1,region_idx{i_region,:}(i_channel),i_unit}];
                            AllRates{i_mouse,i_stim,i_region} = [AllRates{i_mouse,i_stim,i_region};avg_ratesPerBin{1,region_idx{i_region,:}(i_channel),i_unit}];
                        elseif DEC(1,region_idx{i_region,:}(i_channel),i_unit) == 1
                            DecRates{i_mouse,i_stim,i_region} = [DecRates{i_mouse,i_stim,i_region};avg_ratesPerBin{1,region_idx{i_region,:}(i_channel),i_unit}];
                            AllRates{i_mouse,i_stim,i_region} = [AllRates{i_mouse,i_stim,i_region};avg_ratesPerBin{1,region_idx{i_region,:}(i_channel),i_unit}];
                        elseif NC(1,region_idx{i_region,:}(i_channel),i_unit) == 1
                            NCRates{i_mouse,i_stim,i_region} = [NCRates{i_mouse,i_stim,i_region};avg_ratesPerBin{1,region_idx{i_region,:}(i_channel),i_unit}];
                            AllRates{i_mouse,i_stim,i_region} = [AllRates{i_mouse,i_stim,i_region};avg_ratesPerBin{1,region_idx{i_region,:}(i_channel),i_unit}];
                        end
                    end
                end
            end
            
            % Plot the firing rates
            OP.handle = figure(1);
            OP.alpha = 0.2;
            OP.line_width = 2;
            OP.error = 'sem';
            t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
            for i_region = 1:4
                Rates_cell = {IncRates{i_mouse,i_stim,i_region},...
                DecRates{i_mouse,i_stim,i_region},...
                NCRates{i_mouse,i_stim,i_region},...
                AllRates{i_mouse,i_stim,i_region}};
                for i_rates = 1:4
                    subplot(2,2,i_rates)
                    OP.color_area = color_area{i_rates};
                    OP.color_line = color_line{i_rates};
                    plot_areaerrorbar(Rates_cell{i_rates},OP);
                    hold on;
                    ylimits = ylim;
                    ymax = ylimits(2);
                    ylimits = [0 ymax*1.5];
                    line([24 24],ylimits,'Color','k');
                    line([36 36],ylimits,'Color','k');

                    ylim(ylimits);
                    xlim([0 60]);
                    title(names{i_rates})
                    ylabel('Firing Rate (Hz)');
                    xlabel('seconds');
                end
                suptitle(t_bank{i_region})
                saveas(OP.handle,[data_fold 'FiringRates_' stim{i_stim} '_' t_bank{i_region} '.png'],'png');
                clf;
            end
        end
        
    end
end

for i_stim = 1:size(stim,2)
    % Group Average
    OP.handle = figure(1);
    OP.alpha = 0.2;
    OP.line_width = 2;
    OP.error = 'sem';
    t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
    for i_region = 1:4
        Rates_cell = {vertcat(IncRates{:,i_stim,i_region}),...
            vertcat(DecRates{:,i_stim,i_region}),...
            vertcat(NCRates{:,i_stim,i_region}),...
            vertcat(AllRates{:,i_stim,i_region})};
        for i_rates = 1:4
            subplot(2,2,i_rates)
            OP.color_area = color_area{i_rates};
            OP.color_line = color_line{i_rates};
            plot_areaerrorbar(Rates_cell{i_rates},OP);
            hold on;
            ylimits = ylim;
            ymax = ylimits(2);
            ylimits = [0 ymax*1.5];
            line([24 24],ylimits,'Color','k');
            line([36 36],ylimits,'Color','k');

            ylim(ylimits);
            xlim([0 60]);
            title(names{i_rates})
            ylabel('Firing Rate (Hz)');
            xlabel('seconds');
        end
        suptitle(t_bank{i_region})
        saveas(OP.handle,[storage stim{i_stim} slash 'FiringRates_' stim{i_stim} '_' t_bank{i_region} '.png'],'png');
        clf;
    end
end



