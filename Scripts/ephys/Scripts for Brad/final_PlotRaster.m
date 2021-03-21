% Plot the raster plots for different neurons


clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Grad School\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);
Groups = {'Young_VGAT'};
Regions = {'BF','iCT','cCT','Som','ZI'};
Frequencies = [20];

fileTag = '_filt';

%%
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    if ismember(group,Groups)
        if ismember(region,Regions)
            if ismember(mouseList{i,4},Frequencies)
                Filename1 = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz%s.mat',region,frequency,fileTag));
                Filename2 = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_spikeTest.mat',region,frequency));
                if exist(Filename1,'file')
                    fprintf('Loading: %s\n',Filename1);
                    load(Filename1);
                    fprintf('Loading: %s\n',Filename2);
                    load(Filename2);
                    
                    k = 1;
                    
                    % Loop through the different units and generate the plots
                    for i_channel = 1:16     %size(Rates,2)
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
    end
end

saveFilename = ('Raster_INC_070719_C_Som20Hz_ch9k3.eps');
Fig1 = gcf();
print(Fig1,saveFilename,'-depsc','-painters');