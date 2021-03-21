% Look at the inter-event intervals for SLOs

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};
namelist = {'low','high'};


%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
All_Data = cell(length(Groups),length(Regions),2);
for i = 1:size(mouseList,1)
    group = mouseList{i,1};
    mouse = mouseList{i,2};
    region = mouseList{i,3};
    frequency = num2str(mouseList{i,4});
    
    Filename = fullfile(working_dir,mouse,sprintf('BF_%s_%sHz_countSLOs.mat',region,frequency));
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
                
        SLOcount1 = SLOcount(:,1:16);
        SLOcount2 = SLOcount(:,17:32);
        SLOidx1 = SLOidx(:,1:16);
        SLOidx2 = SLOidx(:,17:32);
        
        idx = find(SLOcount1 > 1);
        if ~isempty(idx)
            data = [];
            for i_idx = 1:length(idx)
                current_SLOs = SLOidx1{idx(i_idx)};
                for i_SLO = 2:length(current_SLOs)
                    d = (current_SLOs(i_SLO) - current_SLOs(i_SLO-1))/1000;
                    data = [data,d];
                end
            end
            All_Data{i_group,1,i_frequency} = [All_Data{i_group,1,i_frequency} , data];
        end
        idx = find(SLOcount2 > 1);
        if ~isempty(idx)
            data = [];
            for i_idx = 1:length(idx)
                current_SLOs = SLOidx2{idx(i_idx)};
                for i_SLO = 2:length(current_SLOs)
                    d = (current_SLOs(i_SLO) - current_SLOs(i_SLO-1))/1000;
                    data = [data,d];
                end
            end
            All_Data{i_group,i_region,i_frequency} = [All_Data{i_group,i_region,i_frequency} , data];
        end
    end
end

for i_freq = 1:2
    clf;
    p = 1;
    for i_group = 1:size(All_Data,1)
        fprintf('g: %d\n',i_group);
        group = Groups{i_group};
        group = strsplit(group,'_');
        if group{2}(1) == 'C'
            tmp_freq = '6Hz';
        else
            tmp_freq = '10Hz';
        end
        for i_region = 1:size(All_Data,2)
            fprintf('  r: %d\n',i_region);
            subplot(5,5,p);
            tempData = All_Data{i_group,i_region,i_freq};
            if isempty(tempData)
                set(gca,'Visible','off');
                ax = gca;
                ax.Title.Visible = 'on';
            else
                histogram(tempData,[0:1:20]);%,'Normalization','probability');
                ylimit = ylim;
                xlabel(sprintf('n=%d   avg=%.1f   sem=%.1f ',length(tempData),mean(tempData),std(tempData)/sqrt(length(tempData))));
            end
            if i_region == 1
                ylabel(Groups{i_group},'Interpreter','none');
            end
            if p <= 5
                title(Regions{i_region});
            end
            p = p+1;
        end
    end
    % Save the plot
    SaveFilename = fullfile(working_dir,'Figures_Group',sprintf('PlotSLOs_IEI_%s.png',namelist{i_freq}));
    fprintf('Saving: %s\n',SaveFilename);
    saveas(Fig1,SaveFilename,'png');
%     pause;
end