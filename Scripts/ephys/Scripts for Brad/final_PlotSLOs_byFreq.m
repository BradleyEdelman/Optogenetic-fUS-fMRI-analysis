% Compare the SLO counts at different frequencies within animal

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

[~,~,raw] = xlsread('C:\Users\Lab User\Dropbox\Lee Lab\Data\EphysMouseSummary.xlsx');
rows = 2:137;
mouseList = raw(rows,1:4);

Groups = {'Old_ChAT','Old_ChAT_APP','Young_ChAT','Young_VGAT','Young_VGLUT2'};
Regions = {'BF','iCT','cCT','Som','ZI'};

pThresh = 0.01;

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
All_Data = cell(length(Groups),length(Regions),2,2);
mouseTracker = {};
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
        
        i_mouse = find(cellfun(@(x) isequal([mouse region],x),mouseTracker),1);
        if i == 1
            i_mouse = 1;
        end
        if isempty(i_mouse)
            mouseTracker{end+1} = [mouse region];
            i_mouse = length(mouseTracker);
        end
        All_Data{i_group,1,i_frequency,i_mouse} = SLOcount(:,1:16);
        All_Data{i_group,i_region,i_frequency,i_mouse} = SLOcount(:,17:32);
    end
end

clf;
p = 1;
allINC = 0;
allDEC = 0;
allNC = 0;
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
        nINC = 0;
        nDEC = 0;
        nNC = 0;
        if all(cellfun('isempty',All_Data(i_group,i_region,1,:)))
            set(gca,'Visible','off');
            ax = gca;
            ax.Title.Visible = 'on';
        else
            for i_mouse = 1:size(All_Data,4)
                if ~isempty(All_Data{i_group,i_region,1,i_mouse})
                    fprintf('    %d\n',i_mouse);
                    for i_channel = 1:16
                        if ~isnan(All_Data{i_group,i_region,1,i_mouse}(1,i_channel))
                            lowfreq = All_Data{i_group,i_region,1,i_mouse}(:,i_channel);
                            highfreq = All_Data{i_group,i_region,2,i_mouse}(:,i_channel);
                            
                            lowfreq(isnan(lowfreq)) = [];
                            highfreq(isnan(highfreq)) = [];
                                            
                            pInc = ranksum(lowfreq,highfreq,'Tail','left');
                            pDec = ranksum(lowfreq,highfreq,'Tail','right');
                            
                            if pInc < pThresh
                                lineStyle = '-r';
                                nINC = nINC+1;
                                fprintf('      %d - INC\n',i_channel);
                            elseif pDec < pThresh
                                lineStyle = '-b';
                                nDEC = nDEC+1;
                                fprintf('      %d - DEC\n',i_channel);
                            else
                                lineStyle = '--k';
                                nNC = nNC+1;
                            end
                            
                            n_trials = size(lowfreq,1);
                            both_avg = [mean(lowfreq,1) mean(highfreq,1)];
                            both_std = [std(lowfreq) std(highfreq)];
                            both_sem = [both_std(1)/sqrt(n_trials) both_std(2)/sqrt(n_trials)];
                            
                            plot([0.5 2.5],[0 0],'k','LineWidth',1);
                            hold on;
                            errorbar([1 2],both_avg,both_sem,lineStyle);
                            
                        end
                    end
                end
            end
            
            nTotal = nINC + nDEC + nNC;
            ylimit = ylim;
            text(2.1,0.8*ylimit(2),sprintf('I: %d',nINC));
            text(2.1,0.6*ylimit(2),sprintf('D: %d',nDEC));
            text(2.1,0.4*ylimit(2),sprintf('N: %d',nNC));
            text(2.1,0.2*ylimit(2),sprintf('T: %d',nTotal));
            
            set(gca,'XTick',[1:2]);
            set(gca,'XTickLabel',{tmp_freq,'20Hz'});
            xlim([0.5 2.5]);

        end
        if i_region == 1
            ylabel(Groups{i_group},'Interpreter','none');
        end
        if p <= 5
            title(Regions{i_region});
        end
        p = p+1;
        allINC = allINC + nINC;
        allDEC = allDEC + nDEC;
        allNC = allNC + nNC;
    end
end

% Save the plot
SaveFilename = fullfile(working_dir,'Figures_Group','PlotSLOs_byFreq.png');
fprintf('Saving: %s\n',SaveFilename);
saveas(Fig1,SaveFilename,'png');