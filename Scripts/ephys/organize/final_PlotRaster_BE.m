function final_PlotRaster_BE(storage,base_fold,slash)
% Plot the raster plots for different neurons

%%
stim = {'0_1' '0_5' '1_0'};
figure;

for i_mouse = 1:size(base_fold,1)

    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '_adj.mat'];
        data_file_spikes = [data_fold stim{i_stim} '_spike_test.mat'];
        
        if exist(data_file,'file') && exist(data_file_spikes,'file')
            load(data_file)
            load(data_file_spikes)
                    
            k = 1;
            % Loop through the different units and generate the plots
            for i_channel = 1:size(Rates,2)
                for i_unit = 1:size(Rates,3)
                    if ~isempty(std_ratesPerBin{1,i_channel,i_unit})

                        if (INC(1,i_channel,i_unit) == 1 && (avg_ratesPerPeriod{1,i_channel,i_unit}(2) > 10))
                            fprintf('%d - hit\n',k);

                            subplot(2,1,1);
                            for i_trial=1:size(Spike_TS_adj,1)
                                y = Spike_TS_adj{i_trial,i_channel,i_unit};
                                scatter(y,repmat(i_trial,[1 length(y)]),'+k');
                                ylim([0 21]);
                                yticks([1:20]);
                                set(gca,'xtick',[]);
                                hold on;
                                plot(0:60,repmat(i_trial,61),'-k');
                            end
                            ylabel('Trial #');
                            xlabel('seconds');
                            subplot(2,1,2);

                            % Re-bin the data
                            rebin = 2;
                            temp_data = avg_ratesPerBin{1,i_channel,i_unit};
                            temp_data = reshape(temp_data,rebin,[]);
                            temp_data = mean(temp_data,1);

                            bar(temp_data);
                            hold on;
                            errorbar([10 30 50]/rebin,avg_ratesPerPeriod{1,i_channel,i_unit},sem_ratesPerPeriod{1,i_channel,i_unit},'r');
                            ylabel('Firing Rate(Hz)');
                            xlabel('seconds');
                            xlim([0.5 0.5+60/rebin]);
                            set(gca,'XTick',[0:10:60/rebin]);

                            pause;
                        else
                            fprintf('%d - skip\n',k);
                        end
                        k = k+1;
                        clf;
                    end
                end
            end
            
            
        end
    end
end

saveFilename = ('Raster_INC_070719_C_Som20Hz_ch9k3.eps');
Fig1 = gcf();
print(Fig1,saveFilename,'-depsc','-painters');