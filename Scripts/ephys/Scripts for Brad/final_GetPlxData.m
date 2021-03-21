% Reads in the spike and laser timestamp info as well as the LFP traces
% into matlab from a .plx or .pl2 file
% Using functions from the Matlab Offline Files SDK from Plexon (must be in
% your path)
% Adapted from readall.m in that same SDK

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

% Mice = {'063019_A','063019_B','063019_C','070119_A','070119_B','070119_C','070619_A','070619_B','070619_C','070719_A','070719_C','070719_D','070819_A','070819_B','070819_C'};
Mice = {'072419_A'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'6','10','20'};
Recordings = {'A','B','C','D'};
FileExtension = '.pl2';
overwrite = 0; % 0 to just process new data, 1 for all data

% LFP_TS - Timestamps for LFP
% LFP_mV - Matrix of LFP values in mV for all channels
% Laser_TS - Timestamps for the laser
% Spike_TS - Matrix of spike firing timestamps for all units (Trial [after merging later in pipeline] x Channel x Unit)
% Waveform_mV - Matrix of spike waveforms for all units

%% Perform the operations and save the data
for i_mouse = 1:length(Mice)
    mouse = Mice{i_mouse};
    
    for i_region = 1:length(Regions)
        region = sprintf('BF_%s',Regions{i_region});
        
        for i_freq = 1:length(Frequencies)
            frequency = Frequencies{i_freq};
            
            for i_recording = 1:length(Recordings)
                recording = Recordings{i_recording};
                StartingFileName = fullfile(working_dir,mouse,sprintf('%s_%sHz_%s%s',region,frequency,recording,FileExtension));
                if exist(StartingFileName,'file')
                    SaveFilename = fullfile(working_dir,mouse,sprintf('%s_%sHz_%s.mat',region,frequency,recording));
                    if overwrite == 1 || ~exist(SaveFilename,'file')
                        
                        [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(StartingFileName);
                        
                        disp(['Opened File Name: ' OpenedFileName]);
%                         disp(['Version: ' num2str(Version)]);
                        disp(['Frequency : ' num2str(Freq)]);
                        disp(['Comment : ' Comment]);
                        disp(['Date/Time : ' DateTime]);
                        disp(['Duration : ' num2str(Duration)]);
%                         disp(['Num Pts Per Wave : ' num2str(NPW)]);
%                         disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);
                        % some of the information is only filled if the plx file version is >102
                        if ( Version > 102 )
%                             if ( Trodalness < 2 )
%                                 disp('Data type : Single Electrode');
%                             elseif ( Trodalness == 2 )
%                                 disp('Data type : Stereotrode');
%                             elseif ( Trodalness == 4 )
%                                 disp('Data type : Tetrode');
%                             else
%                                 disp('Data type : Unknown');
%                             end
%                             
%                             disp(['Spike Peak Voltage (mV) : ' num2str(SpikePeakV)]);
%                             disp(['Spike A/D Resolution (bits) : ' num2str(SpikeADResBits)]);
%                             disp(['Slow A/D Peak Voltage (mV) : ' num2str(SlowPeakV)]);
%                             disp(['Slow A/D Resolution (bits) : ' num2str(SlowADResBits)]);
                        end
                        
                        
                        % Get LFP Data
                        % This is complicated by channel numbering.
                        % The number of samples for analog channel 0 is stored at slowcounts(1).
                        % Note that analog ch numbering starts at 0, not 1 in the data, but the
                        [nad,adfreqs] = plx_adchan_freqs(OpenedFileName);
                        [nad,adgains] = plx_adchan_gains(OpenedFileName);
                        [nad,adnames] = plx_adchan_names(OpenedFileName);
                        nLFP = sum(adnames(:,1) == 'F'); % F for field potential, to get just LFP
                        % preallocate for speed
                        LFP_mV = cell(1,nLFP);
                        k = 1;
                        for ich = find(adnames(:,1) == 'F')'
                            [adfreq, nad, tsad, fnad, LFP_mV{k}] = plx_ad_v(OpenedFileName, ich-1);
                            k = k + 1;
                        end
                        if length(tsad) > 1
                            warning('This file contains %d fragments in the LFP timestamps (expected 1)',length(tsad));
                        end
                        LFP_TS = [];
                        for i_tsad = 1:length(tsad)
                            LFP_TS = [LFP_TS tsad(i_tsad):1/1000:((fnad(i_tsad)-1)/1000+tsad(i_tsad))];
                        end
                        
                        % Get Laser Timing Data
                        % need the event chanmap to make any sense of these
                        [nev,evnames] = plx_event_names(OpenedFileName);[u,evchans] = plx_event_chanmap(OpenedFileName);
                        evch = evchans(1);
                        if (evch == 257)
                            [nevs, Laser_TS, svStrobed] = plx_event_ts(OpenedFileName, evch);
                        else
                            [nevs, Laser_TS, svdummy] = plx_event_ts(OpenedFileName, evch);
                        end
                        
                        
                        % get some info about the spike channels
                        [nspk,spk_filters] = plx_chan_filters(OpenedFileName);
                        [nspk,spk_gains] = plx_chan_gains(OpenedFileName);
                        [nspk,spk_threshs] = plx_chan_thresholds(OpenedFileName);
                        [nspk,spk_names] = plx_chan_names(OpenedFileName);
                        [tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);
                        % tscounts, wfcounts are indexed by (unit+1,channel+1)
                        % tscounts(:,ch+1) is the per-unit counts for channel ch
                        % sum( tscounts(:,ch+1) ) is the total wfs for channel ch (all units)
                        % [nunits, nchannels] = size( tscounts )
                        % To get number of nonzero units/channels, use nnz() function
                        % gives actual number of units (including unsorted) and actual number of
                        % channels plus 1
                        [nunits1, nchannels1] = size(tscounts);
                        nunits = nunits1 - 1;
                        nchannels = nchannels1 - 1;
                        % preallocate for speed
                        Spike_TS = cell(nunits, nchannels);
                        Waveform_mV = cell(nunits, nchannels);
                        for iunit = 1:nunits   % sorted units start at 1, 0 is for unsorted units
                            for ich = 1:nchannels
                                if (tscounts(iunit, ich+1) > 0)
                                    % get the timestamps and waveforms for this channel and unit
                                    [wave_n, wave_npw, Spike_TS{iunit,ich}, Waveform_mV{iunit,ich}] = plx_waves_v(OpenedFileName, ich, iunit-1);
                                end
                            end
                        end
                        Spike_TS = permute(Spike_TS,[3 2 1]);
                        Waveform_mV = permute(Waveform_mV,[3 2 1]);
                        
                        % Save the Data to a .mat file
                        fprintf('Saving...\n\n');
                        save(SaveFilename,'LFP_TS','LFP_mV','Spike_TS','Waveform_mV','Laser_TS');
                    end
                end
            end
        end
    end
end