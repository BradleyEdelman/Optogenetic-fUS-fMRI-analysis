% Plot the delays between the laser onset and spikes

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

xlimit = [-0.001 0.030];
hist_binsize = 0.000025; %seconds
namelist = {'INC_low','INC_high','DEC_low','DEC_high','NC_low','NC_high','ALL_low','ALL_high','INC_both','DEC_both','NC_both','ALL_both'};
coord = {[1,1],[1,2],[2,1],[2,2],[3,1],[3,2],[4,1],[4,2],[1,3],[2,3],[3,3],[4,3]}; % For indexing results for namelist

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
hold on;
for curr_ch = 1:32
    All_Data = cell(length(Groups),length(Regions),3,4);
    for i = 1:size(mouseList,1)
        group = mouseList{i,1};
        mouse = mouseList{i,2};
        region = mouseList{i,3};
        frequency = num2str(mouseList{i,4});

        Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_SpikeDelays.mat',region,frequency));
        Filename2 = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_spikeTest.mat',region,frequency));
        if exist(Filename,'file')
            fprintf('Loading: %s\n',Filename);
            load(Filename);
            load(Filename2);

            [~,i_group] = ismember(group,Groups);
            [~,i_region] = ismember(region,Regions);
            if frequency == '20'
                i_frequency = 2;
            else
                i_frequency = 1;
            end

            % Calculate values for the spike delays
            inc = cell(size(Spike_Delays,1),32,size(Spike_Delays,3));
            dec = inc;
            nc = inc;
            all = inc;
            for i_ch = curr_ch
                for i_unit = 1:size(Spike_Delays,3)
                    if INC(1,i_ch,i_unit) == 1
                        inc(:,i_ch,i_unit) = Spike_Delays(:,i_ch,i_unit);
                    elseif DEC(1,i_ch,i_unit) == 1
                        dec(:,i_ch,i_unit) = Spike_Delays(:,i_ch,i_unit);
                    elseif NC(1,i_ch,i_unit) == 1
                        nc(:,i_ch,i_unit) = Spike_Delays(:,i_ch,i_unit);
                    end
                    all(:,i_ch,i_unit) = Spike_Delays(:,i_ch,i_unit);
                end
            end        
            Data = {inc,dec,nc,all};
            if curr_ch <= 16
                i_region = 1;
            end
            for i_flat = 1:length(Data)
                input = Data{i_flat};
                input = reshape(permute(input,[1 3 2]),[size(input,1)*size(input,3),size(input,2)]);

                All_Data{i_group,i_region,i_frequency,i_flat} = [All_Data{i_group,i_region,i_frequency,i_flat};cell2mat(reshape(input(:,:),[],1))];
                All_Data{i_group,i_region,3,i_flat} = [All_Data{i_group,i_region,3,i_flat};cell2mat(reshape(input(:,:),[],1))];
            end
            clear Data;
        end
    end

    for ii = 1:length(namelist)
        p = 1;
        for i_group = 1:size(All_Data,1)
            for i_region = 1:size(All_Data,2)
                subplot(5,5,p);
                if ~isempty(All_Data{i_group,i_region,coord{ii}(2),coord{ii}(1)})

                    % Plot the data
                    pts = All_Data{i_group,i_region,coord{ii}(2),coord{ii}(1)}; 

                    pts = pts(pts >= 1/40000); % Get rid of values less than the sampling frequency

                    histogram(pts,[0:hist_binsize:0.200]);%,'Normalization','probability');
                    hold on;
                    ylimit = ylim;
                    cla;
                    verts = [.005,.01,.015,.02,.025];
                    for i_v = 1:length(verts)
                        p1 = line([verts(i_v) verts(i_v)],[ylimit(1) ylimit(2)],'Color','k');
                        p1.Color(4) = 0.3;
                    end
                    histogram(pts,[0:hist_binsize:0.200]);%,'Normalization','probability');


                    xticks = [0:0.005:0.20];
                    set(gca,'XTick',xticks);
                    set(gca,'XTickLabel',xticks*1000);
                    xlim(xlimit);

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
        % Save the plot
        fprintf('Saving: %s',fullfile(working_dir,'Figures_Group','PlotSpikeDelays_byChannel',sprintf('SpikeDelays_Ch%d_%s.png\n',curr_ch,namelist{ii})));
        saveas(Fig1,fullfile(working_dir,'Figures_Group','PlotSpikeDelays_byChannel',sprintf('SpikeDelays_Ch%d_%s.png',curr_ch,namelist{ii})),'png');
%         pause;
        clf;
    end
end