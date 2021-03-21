function final_PlotStackedBar_BE(storage,base_fold,slash)
% Plot the stacked bar graph for the different neuronal responses

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');

stim = {'0_1' '0_5' '1_0'};

hold on;
% Data = cell(length(Groups),length(Regions),2);
for i_mouse = 1:size(base_fold,1)
    
    for i_stim = 1:size(stim,2)
        data_fold = [storage base_fold{i_mouse} slash stim{i_stim} slash];
        data_file = [data_fold stim{i_stim} '_spike_test.mat'];
        
        if exist(data_file,'file')
            load(data_file);
            
            % Calculate values for the stacked bar graphs for each region
            % at the group level
            E1_nINC = sum(sum(INC(:,1:8,:)));
            E2_nINC = sum(sum(INC(:,9:16,:)));
            E3_nINC = sum(sum(INC(:,17:24,:)));
            E4_nINC = sum(sum(INC(:,25:32,:)));
            
            E1_nDEC = sum(sum(DEC(:,1:8,:)));
            E2_nDEC = sum(sum(DEC(:,9:16,:)));
            E3_nDEC = sum(sum(DEC(:,17:24,:)));
            E4_nDEC = sum(sum(DEC(:,25:32,:)));
            
            E1_nNC = sum(sum(NC(:,1:8,:)));
            E2_nNC = sum(sum(NC(:,9:16,:)));
            E3_nNC = sum(sum(NC(:,17:24,:)));
            E4_nNC = sum(sum(NC(:,25:32,:)));
            
            E1_total = E1_nINC + E1_nDEC + E1_nNC;
            E2_total = E2_nINC + E2_nDEC + E2_nNC;
            E3_total = E3_nINC + E3_nDEC + E3_nNC;
            E4_total = E4_nINC + E4_nDEC + E4_nNC;
            
            Data{i_mouse,i_stim,1} = [E1_nINC, E1_nDEC, E1_nNC];
            Data{i_mouse,i_stim,2} = [E2_nINC, E2_nDEC, E2_nNC];
            Data{i_mouse,i_stim,3} = [E3_nINC, E3_nDEC, E3_nNC];
            Data{i_mouse,i_stim,4} = [E4_nINC, E4_nDEC, E4_nNC];
        end
    end
end

t_bank = {'L CPu' 'L Cortex' 'R CPu' 'R Cortex'};
for i_region = 1:4
    
    for i_stim = 1:size(stim,2)
        if ~isempty(vertcat(Data{:,i_stim,i_region})) 
            
            total = vertcat(Data{:,i_stim,i_region});
            pctInc(i_stim,i_region) = sum(total(:,1))/sum(total(:))*100;
            pctDec(i_stim,i_region) = sum(total(:,2))/sum(total(:))*100;
            pctNC(i_stim,i_region) = sum(total(:,3))/sum(total(:))*100;
        else
            pctInc(i_stim,i_region) = 0;
            pctDec(i_stim,i_region) = 0;
            pctNC(i_stim,i_region) = 0;
        end
    end
    
    % Create the stacked bar graphs
    figure(Fig1)
    subplot(2,2,i_region)
    BAR = [pctInc(1,i_region),pctDec(1,i_region),pctNC(1,i_region);...
        pctInc(2,i_region),pctDec(2,i_region),pctNC(2,i_region);...
        pctInc(3,i_region),pctDec(3,i_region),pctNC(3,i_region)];
    B = bar(BAR,'stacked');
    B(1).FaceColor = 'r';
    B(2).FaceColor = 'b';
    B(3).FaceColor = [128 128 128]/255;
    set(gca,'XTick',[1 2 3]);
    set(gca,'TickDir','out'); box off;
    set(gca,'XTickLabel',{'0.1' '0.5' '1.0'});
    xlim([0.5 3.5]);
    set(gca,'YTick',[0:20:100]);
    ylim([0 100]);
    set(gca,'color','none');
    title(t_bank{i_region})
end
% Save the plot
saveas(Fig1,[storage 'Spike_Change_Stacked_Bars_Grp.png'],'png');