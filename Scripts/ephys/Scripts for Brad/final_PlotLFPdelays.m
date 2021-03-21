% Plot the delays between the laser onset and first LFP peak

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
Data = cell(length(Groups),length(Regions),2);
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_LFPdelays.mat',region,frequency));
    if exist(Filename,'file')
        fprintf('Loading: %s\n',Filename);
        load(Filename);
        
        % Calculate values for the LFP delays
        Delays = mean(LFP_Delays,1);
        if size(Delays,2) < 32
            new = zeros(1,32,size(Delays,3));
            for i_ch = 1:size(Delays,2)
                new(:,i_ch,:) = Delays(:,i_ch,:);
            end
            Delays = new;
        end
        
        [~,i_group] = ismember(group,Groups);
        [~,i_region] = ismember(region,Regions);
        if frequency == '20'
            i_frequency = 2;
        else
            i_frequency = 1;
        end
        Data{i_group,1,i_frequency} = [Data{i_group,1,i_frequency};Delays(:,1:16,:)];
        Data{i_group,i_region,i_frequency} = [Data{i_group,i_region,i_frequency};Delays(:,17:end,:)];
    end
end

for i_thresh = 1:length(maxThresh)
    p = 1;
    for i_group = 1:size(Data,1)
        for i_region = 1:size(Data,2)
            subplot(5,5,p);
            if ~isempty(Data{i_group,i_region,1})
                
                % Plot the data
                low_freq = Data{i_group,i_region,1}(:,:,i_thresh);
                high_freq = Data{i_group,i_region,2}(:,:,i_thresh);
                
                scatter(reshape(repmat([1:16],[size(low_freq,1),1]),[],1),reshape(low_freq,[],1),'xr');
                hold on;
                scatter(reshape(repmat([1:16],[size(high_freq,1),1]),[],1),reshape(high_freq,[],1),'xb');
                
                set(gca,'XTick',[1:16]);
                set(gca,'XTickLabel',[17:32],'fontsize',8);
                xlim([0 17]);
%                 ylim([0 200]);
                if i_region == 1
                    ylabel(Groups{i_group},'Interpreter','none');
                    set(gca,'XTickLabel',[1:16],'fontsize',8);
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
    fprintf('Saving: %s\n',fullfile(working_dir,'Figures_Group',sprintf('LFPdelays_%dpct.png',maxThresh(i_thresh))));
    saveas(Fig1,fullfile(working_dir,'Figures_Group',sprintf('LFPdelays_%dpct.png',maxThresh(i_thresh))),'png');
%     pause;
    clf;
end