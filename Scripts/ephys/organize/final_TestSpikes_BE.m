function final_TestSpikes_BE(storage,base_fold,slash)
% Test the spiking data for significant changes, save an output file with
% the results

binsize = 1.00; % seconds
pThresh = 0.05;
ttest_flag = 0; % If 1, perform t-test instead of sign-rank test

%%
% 24 sec pre, 12 sec peri, 24 sec post 
preIdx = 1:24;
stimIdx = 24:36;
postIdx = 47:60;

Stim = {'0_1' '0_5' '1_0'};
for i_mouse = 1:size(base_fold,1)
   
    for i_stim = 1:size(Stim,2)
        data_fold = [storage base_fold{i_mouse} slash Stim{i_stim} slash];
        data_file = [data_fold Stim{i_stim} '_adj.mat'];
        save_file = [data_fold Stim{i_stim} '_spike_test.mat'];
            
        if exist(data_file,'file')
            load(data_file);
        
            n_trials = size(Spike_TS_adj,1);
            n_channels = size(Spike_TS_adj,2);
            n_units = size(Spike_TS_adj,3);
        
            Rates = cell(n_trials,n_channels,n_units);
            avg_ratesPerBin = cell(1,n_channels,n_units);
            std_ratesPerBin = cell(1,n_channels,n_units);
            sem_ratesPerBin = cell(1,n_channels,n_units);
            avg_ratesPerPeriod = cell(1,n_channels,n_units);
            std_ratesPerPeriod = cell(1,n_channels,n_units);
            sem_ratesPerPeriod = cell(1,n_channels,n_units);
        
            INC = zeros(1,n_channels,n_units); % 1 where there is an increase
            DEC = zeros(1,n_channels,n_units); % 1 where there is a decrease
            NC = zeros(1,n_channels,n_units); % 1 where there is no change
        
        
            % Bin the data
            fprintf('Binning & Testing Spike Data, Bin Size = %dms\n',binsize*1000);
            for i_channel = 1:n_channels
                for i_unit = 1:n_units

                    continue_flag = 0;
                    pre = zeros(n_trials,length(preIdx));
                    stim = zeros(n_trials,length(stimIdx));
                    post = zeros(n_trials,length(postIdx));

                    for i_trial = 1:n_trials
                        if ~isempty(Spike_TS_adj{i_trial,i_channel,i_unit})
                            binned = histcounts(Spike_TS_adj{i_trial,i_channel,i_unit},[0:binsize:60]);
                            Rates{i_trial,i_channel,i_unit} = binned;
                            if continue_flag == 0
                                continue_flag = 1;
                            end
                        end
                    end

                    if continue_flag == 1
                        % Calculate the mean & variance across trials within each bin
                        ratesPerBin = cell2mat(Rates(:,i_channel,i_unit));
                        avg_ratesPerBin{1,i_channel,i_unit} = mean(ratesPerBin,1);
                        std_ratesPerBin{1,i_channel,i_unit} = std(ratesPerBin);
                        sem_ratesPerBin{1,i_channel,i_unit} = std(ratesPerBin)/sqrt(n_trials);

                        % Calculate the mean & variance across trials within each period (pre, stim, post)
                        pre = mean(ratesPerBin(:,preIdx),2);
                        stim = mean(ratesPerBin(:,stimIdx),2);
                        post = mean(ratesPerBin(:,postIdx),2);
                        avg_ratesPerPeriod{1,i_channel,i_unit} = [mean(pre,1) mean(stim,1) mean(post,1)];
                        std_ratesPerPeriod{1,i_channel,i_unit} = [std(pre) std(stim) std(post)];
                        sem_ratesPerPeriod{1,i_channel,i_unit} = [std(pre) std(stim) std(post)]/sqrt(n_trials);

                        % Test the data
                        if ttest_flag == 1
                            [~,pInc] = ttest(pre,stim,'Tail','left');
                            [~,pDec] = ttest(pre,stim,'Tail','right');
                        else
                            pInc = signrank(pre,stim,'Tail','left');
                            pDec = signrank(pre,stim,'Tail','right');
                        end
                        
                        if pInc < pThresh
                            INC(1,i_channel,i_unit) = 1; %q = (pInc < 0.05) + (pInc < 0.005) + (pInc < 0.001);
                            fprintf('%d %d pause inc\n',i_channel,i_unit);%pause;
                        elseif pDec < pThresh
                            DEC(1,i_channel,i_unit) = 1; %q = (pDec < 0.05) + (pDec < 0.005) + (pDec < 0.001);
                            fprintf('%d %d pause dec\n',i_channel,i_unit);%pause;
                        else
                            NC(1,i_channel,i_unit) = 1; %q = 0;
                        end
                    end
                end
            end
        
        % Save the results
        fprintf('   [Total Inc Dec NC] = %d - %d %d %d\n',sum(sum(INC))+sum(sum(DEC))+sum(sum(NC)),sum(sum(INC)),sum(sum(DEC)),sum(sum(NC)));
        fprintf(['\nSaving: ' save_file '\n']);
        save(save_file,'pThresh','binsize','preIdx','stimIdx','postIdx','Rates','INC','DEC','NC','avg_ratesPerBin','std_ratesPerBin','sem_ratesPerBin','avg_ratesPerPeriod','std_ratesPerPeriod','sem_ratesPerPeriod');
        end
    end
end