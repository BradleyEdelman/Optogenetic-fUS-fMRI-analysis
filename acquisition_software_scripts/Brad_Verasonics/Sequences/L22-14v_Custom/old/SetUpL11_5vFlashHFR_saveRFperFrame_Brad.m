% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: RunSetUpL11_5vFlashHFR_saveRFperFrame.m
% "High Frame Rate" acquisition for ultrafast imaging with saving RF data per frame in realtime.
% This script shows how to use external function to save the RF data of each frame in
% realtime. Here, only non-zeros channel data will be saved to reduce the
% saving time. The default frame limit is set to 100 "superframes". Please change
% to a bigger number for more superframes.
%
% Description:
% To support higher acquisition frame rates with reduced DMA overhead, this script acquires
%   a large set of T/R acquisitions into each RcvBuffer 'super' frame, performing a transferToHost only after
%   each group of "numAcqs" acquisitions.
%   Even though each acquisition is intended to be reconstructed and processed into a display frame, live image reconstruction
%   only processes the first acquisition in the super frame in order not to slow down real-time acquisition,
%   and yet provide visual feedback to help the user collect desired data. To post-process all acquisitions and superframes
%   present in the RcvBuffer, the data is saved in a matfile when VSX quits (see end of the script).
%   Unfortunately, the current system software does not permit running this same script to process all of the acquisitions,
%   even in playback (simulationMode=2), because a sequence can only reconstruct a maximum of 16 acquisitions
%   in a given receive frame, and therefore cannot include all of the data transferred in one DMA (here called a super frame).
%
% To illustrate the HFR acquisition and display in simulation, the 'movePoints' function is invoked for every acquisition,
%   resulting in a very discontinous apparent motion of the scatterers in "real-time".
%
% For convenience, this script is currently set to launch VSX
% automatically. In order to save each frame in a correct order, a "sync"
% command is required for the hardware to wait for the software to finish
% saving RF data, image reconstruction and image display.
% Therefore, a warning message "timeToNextAcq duration too short" might occur
% if the interval (SeqControl(5).argument in this script) between two frames is
% not long enough.
%
% Last update:
% 05/23/2016 - test with SW 3.0.7

clear all

% --- Frequently Modified Parameters ------------------------------------------
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 128;   % This should preferrably be a multiple of 128 samples.

P.numAcqs = 200;      % no. of Acquisitions in a Receive frame (this is a "superframe")
simulateMode = 0;   % set to acquire data using Vantage 128 hardware

P.path = 'C:\Users\verasonics\Documents\Vantage-3.4.2-1806191000\';
P.filePrefix = 'SaveIQData';

P.time = clock;
P.dateStr = strcat('_',num2str(P.time(2)), '-', num2str(P.time(3)), '-',...
    num2str(P.time(1)));

P.saveAcquisition = 1; %Default doesn't save
P.settingsNumber = 1; %Which version of the settings are you on?
P.settingsChanged = 1; %Whether the settings have changed since the last save.  Starts at 1 so that it doesn't automatically iterate to 2

P.runNumber = 1; %What run on the current setting?
P.itNumber = 1; %What iteration on the current run?
% -----------------------------------------------------------------------------
na = 5;      % Set na = number of angles.
if (na > 1)
    dtheta = (12*pi/180)/(na-1);
    P.startAngle = -12*pi/180/2;
else
    dtheta = 0;
    P.startAngle=0;
end

PRF=na*500; % Acquisition Rate

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;

% Specify Trans structure array.
Trans.name = 'L22-14v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 25;      % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [0.4, 0, 0.25];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 1024*P.numAcqs*na;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 1;       % number of 'super frames'

Resources.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % only one intermediate buffer needed.
Resource.InterBuffer(1).pagesPerFrame = P.numAcqs;

Resource.ImageBuffer(1).numFrames = 1;

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,4,1];


% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end
for n = 1:na   % na transmit events
    TX(n).Steer = [(P.startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [330 560 780 1010 1023 1023 1023 1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 1), 1, P.numAcqs*na);  % movepoints EVERY acquisition to illustrate superframe concept
                                                                    % real-time images will look "jerky" but using the reconstructAll script,
                                                                    % playback process all acquisitions and shows smooth displacement

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numAcqs*na*(i-1) + 1).callMediaFunc = 1;  % move points only once per frame
    for j = 1:P.numAcqs
        % -- Acquisitions for 'super' frame.
        rcvNum = P.numAcqs*na*(i-1) + j;
        Receive(rcvNum).Apod(:)=1;
        Receive(rcvNum).framenum = i;
        Receive(rcvNum).acqNum = j;
    end
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:P.numAcqs*na);

% Define ReconInfo structures.  (just one ReconInfo to process one acquisition per superframe)
ReconInfo = repmat(struct('mode', 'accumIQ', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...  % use the first acquisition of each frame
                   'regionnum', 1), 1, P.numAcqs*na);
               
% % - Set specific ReconInfo attributes.
if na>1
    
    for j = 1:P.numAcqs*na  % For each row in the column
        
        if isequal(rem(j,na),1)
            ReconInfo(j).mode = 'replaceIQ'; % replace IQ data every na acquisitions
        end
        
        if isequal(rem(j,na),0) % cycle through angle acquisitions in each page
            ReconInfo(j).txnum=na;
        else
            ReconInfo(j).txnum=rem(j,na);
        end
        
        ReconInfo(j).rcvnum = j;
        ReconInfo(j).pagenum = ceil(j/na);
        
        if isequal(rem(j,na),0) % Acummulate angle data every na acquisitions
            ReconInfo(j).mode = 'accumIQ_replaceIntensity';
        end
        
    end
        
else
    ReconInfo(1).mode = 'replaceIntensity';
end               

% Save Data Process
Process(1).classname = 'External';
Process(1).method = 'saveIQData'; % calls the 'saveRFperFrame' function
Process(1).Parameters = {'srcbuffer','inter',... % name of buffer to process.
        'srcbufnum',1,...
        'srcframenum',-1,... % process the most recent frame.
        'dstbuffer','none'};
    
Process(2).classname = 'External';
Process(2).method = 'saveRFperFrame'; % calls the 'saveRFperFrame' function
Process(2).Parameters = {'srcbuffer','receive',... % name of buffer to process.
        'srcbufnum',1,...
        'srcframenum',-1,... % process the most recent frame.
        'dstbuffer','none'};    
    
    
% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between acquisitions
SeqControl(2).argument = 400;  % 200 usecs
SeqControl(3).command = 'timeToNextAcq';  % time between superframes
SeqControl(3).argument = 400;  % 20 ms
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'sync';
SeqControl(5).argument = 0.5e6;   % .5s
nsc = 6; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numAcqs*na
        Event(n).info = 'Acquire RF';
        
        if isequal(rem(j,na),0)
            Event(n).tx = na;
        else
            Event(n).tx=rem(j,na);
        end
        
        Event(n).rcv = P.numAcqs*na*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [3,nsc];
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        nsc = nsc + 1;

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 5;
    n = n+1;

    % Do reconstruction and processing for 1st sub frame
%     Event(n).info = 'Reconstruct';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 1;
%     Event(n).process = 0;
%     Event(n).seqControl = 4;
%     n = n+1;
    
%     Event(n).info = 'Save IQ data';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 1;
%     Event(n).seqControl = 0;
%     n=n+1;

end

Event(n).info = 'Save RF data';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 2;
Event(n).seqControl = 0;
n=n+1;

% --- If this last event is not included, the sequence stops after one pass, and enters "freeze" state
%     Pressing the freeze button runs the "one-shot" sequence one more time
%     For live acquisition in mode 0, simply comment out the 'if/end' statements and manually freeze and exit when the data looks good.
if simulateMode==2 || simulateMode==0 %  In live acquisiton or playback mode, run continuously, but run only once for all frames in simulation
    Event(n).info = 'Jump back to first event';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 1;
end

% 
EF(1).Function = text2cell('%EF#1');
EF(2).Function = text2cell('%saveRFperFrame');

% % Specify factor for converting sequenceRate to frameRate.
frameRateFactor = P.numAcqs;

% Save all the structures to a .mat file.
% and invoke VSX automatically
filename = 'L11-5vFlashHFR_saveRFperFrame'; % used to launch VSX automatically
save(['MatFiles/',filename]);
VSX

return


% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%EF#1
saveIQData(IQData)
    
if evalin('base','P.saveAcquisition')
            % File Naming and IQ data
    %Read in the misc variables struct 
    P = evalin('base','P');

    %Now we want to handle the settings file, to make sure we are saving
    %correctly

    %If the settings have changed since the last time, reset the boolean to
    %false so that new changes will propogate.  Also reset the run number
    %since it's the first run on the new settings
    if P.itNumber == 1
        %Check for previous settings
        while exist(strcat(P.path,P.filePrefix,P.dateStr,...
        '_Run-',int2str(P.runNumber),'_It-',int2str(P.itNumber),'_IQ.mat'),'file')
            P.runNumber = P.runNumber+1;
        end

       %TODO: Line to invoke save preSet here.
    end

    %Calculate the file name for any iteration specific file
    fileName = strcat(P.path,P.filePrefix,P.dateStr,...
        '_Run-',int2str(P.runNumber),'_It-',int2str(P.itNumber));

    %Save the IQ data for the run.
    tic
    save(strcat(fileName,'_IQ'),'IQData'); %Save the IQ data     
    fprintf('saving time for frame %g = %g s \n',P.itNumber, toc);
    %Modify the iteration number
    P.itNumber = P.itNumber+1;
    assignin('base','P',P);
end
return
%EF#1


%saveRFperFrame - save RF
saveRFperFrame(RData)

persistent frameNum RFfilename

% file size can be reduced by elimating all zeros
TXApod = evalin('base','TX.Apod');
endSample = evalin('base','Receive(end).endSample');

frameNumLimit = 100;
numLength = max(ceil(log10(abs(frameNumLimit))),1)+1;

if isempty(frameNum)
    frameNum = 1;
    RFfilename = ['RFdata_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')];
else
    frameNum = frameNum + 1;
    if frameNum > frameNumLimit, frameNum = 1; end % set a limit for testing
end

fname = ['frame_',num2str(frameNum,['%0',num2str(numLength),'.0f'])];
newRData = RData(1:endSample,find(TXApod));
eval([fname,' = newRData;']);

tic
if isequal(frameNum,1)
    save(RFfilename,fname,'-v6');
else
    save(RFfilename,fname,'-v6','-append');
end
fprintf('saving time for frame %g = %g s \n',frameNum, toc);
return
%saveRFperFrame