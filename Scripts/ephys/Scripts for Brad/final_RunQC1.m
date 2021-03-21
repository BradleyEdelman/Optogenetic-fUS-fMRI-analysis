% Perform quality control step 1 on the individual recordings

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

Mice = {'072419_A'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'6','10','20'};
Recordings = {'A','B','C','D'};

%%
for i_mouse = 1:length(Mice)
    mouse = Mice{i_mouse};
    for i_region = 1:length(Regions)
        region = sprintf('BF_%s',Regions{i_region});
        for i_freq = 1:length(Frequencies)
            frequency = Frequencies{i_freq};
            QC_Status = [];
            save_flag = 0;
            for i_recording = 1:length(Recordings)
                recording = Recordings{i_recording};
                Filename = fullfile(working_dir,mouse,sprintf('%s_%sHz_%s.mat',region,frequency,recording));
                if exist(Filename,'file')
                    save_flag = 1;
                    msg = sprintf('Loading: %s\n',Filename);
                    fprintf('Loading: %s\n',Filename);
                    load(Filename);
                    QC_Status = [QC_Status msg];
                    
                    % Check for duplicate LFP channels
                    msg = sprintf('   Checking LFP for duplicates\n');
                    fprintf(msg);
                    QC_Status = [QC_Status msg];
                    n_channels = length(LFP_mV);
                    sum_mV = zeros(1,n_channels);
                    for i_channel = 1:n_channels
                        sum_mV(i_channel) = sum(LFP_mV{i_channel});
                    end
                    [~,ia,~] = unique(sum_mV,'stable');
                    if length(ia) == n_channels
                        msg = sprintf('      PASS\n');
                        fprintf(msg);
                        QC_Status = [QC_Status msg];
                    else
                        msg1 = sprintf('      Possible LFP Duplicates Detected on Channels: ');
                        msg2 = sprintf('%d  ',setdiff(1:n_channels,ia));
                        fprintf(msg1);
                        disp(msg2);
                        QC_Status = [QC_Status msg1 msg2 newline];
                    end
                    
                    % Check for duplicate spike units
                    msg = sprintf('   Checking spike units for duplicates\n');
                    fprintf(msg);
                    QC_Status = [QC_Status msg];
                    count_units = 0;
                    drop_0_flag = 0;
                    n_channels = size(Spike_TS,2);
                    n_units = size(Spike_TS,3);
                    sum_TS = zeros(n_channels,n_units);
                    for i_channel = 1:size(Spike_TS,2)
                        for i_unit = 1:size(Spike_TS,3)
                            if ~isempty(Spike_TS{1,i_channel,i_unit})
                                sum_TS(i_channel,i_unit) = sum(Spike_TS{1,i_channel,i_unit});
                                count_units = count_units + 1;
                            elseif drop_0_flag == 0
                                drop_0_flag = 1;
                            end
                        end
                    end
                    [~,ia,~] = unique(sum_TS,'stable');
                    [r,c] = ind2sub([n_channels,n_units],ia);
                    unique_units = length(r);
                    if drop_0_flag == 1
                        unique_units = unique_units - 1;
                    end
                    if unique_units == count_units
                        msg = sprintf('      PASS\n');
                        fprintf(msg);
                        QC_Status = [QC_Status msg];
                    else
                        msg1 = sprintf('      Possible Spike Unit Duplicates Detected on Channels: ');
                        msg2 = sprintf('%d  ',setdiff(1:n_channels,r));
                        fprintf(msg1);
                        disp(msg2);
                        QC_Status = [QC_Status msg1 msg2 newline];
                    end
                    msg = fprintf('\n\n');
                    QC_Status = [QC_Status newline newline];
                end
            end
            
            % Save the results
            if save_flag == 1
                OutputFilename = fullfile(working_dir,mouse,sprintf('QC1_%s_%sHz.txt',region,frequency));
                dlmwrite(OutputFilename,QC_Status,'delimiter','');
            end
        end
    end
end