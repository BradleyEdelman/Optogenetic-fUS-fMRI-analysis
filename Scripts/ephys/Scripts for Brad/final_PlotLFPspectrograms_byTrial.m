% Plot the LFP Spectrogram for individual trials & sites (modified to save individual plots)

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

Mice = {'072419_A'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'6','10','20'};

t_window = 1000;
t_overlap = 800;
nfft = 1000;
fs = 1000;
climit = [-150 -35];
ylimit = [0 80];
i_png = 2641;

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
for i_mouse = 1:length(Mice)
    mouse = Mice{i_mouse};
    for i_region = 1:length(Regions)
        region = sprintf('BF_%s',Regions{i_region});
        for i_freq = 1:length(Frequencies)
            frequency = Frequencies{i_freq};
            Filename = fullfile(working_dir,mouse,sprintf('%s_%sHz_adj.mat',region,frequency));
            if exist(Filename,'file')
                fprintf('Loading: %s\n',Filename);
                load(Filename);
                
                n_channels = size(LFP_mV_adj,2);
                n_trials = size(LFP_mV_adj,1);
                
                % Plot the indivudal spectrograms
                for i_trial = 1:n_trials
                    fprintf('Plotting/Saving: %s  %s  %sHz  Ch. 1-32 Trial %d\n',mouse,region(4:end),frequency,i_trial);
                    for i_channel = 1:n_channels
                        ts = LFP_mV_adj{i_trial,i_channel};
                        subplot(4,8,i_channel);
%                         spectrogram(ts,t_window,t_overlap,nfft,fs,'yaxis');
                        
                        [s,f,t,p] = spectrogram(ts,t_window,t_overlap,nfft,fs,'yaxis','power');
                        if isequal(i_trial,1)
                            imagesc(mag2db(p)); axis xy
                        else
                            g = get(gca,'Children');
                            g.CData = mag2db(p);
                        end
                        
                        ax=gca;
                        ax.CLim = climit;
                        ylim(ylimit);
                        colorbar('off');xlabel('');ylabel('');
                        if i_channel == 1
                            xlabel(sprintf('%s %s %sHz Trial %d',mouse,region(4:end),frequency,i_trial),'Interpreter','none');
                        end
%                         title({sprintf('%s',mouse),sprintf('%s_%sHz',region,frequency),sprintf('Channel %d Trial %d',i_channel,i_trial)},'Interpreter','none');
                    end
                    saveas(Fig1,fullfile(working_dir,'Figures_Indiv','PlotLFPspectrograms_byTrial',sprintf('%d.png',i_png)),'png');
%                     clf;
                    i_png = i_png+1;
                end
            end
        end
    end
end