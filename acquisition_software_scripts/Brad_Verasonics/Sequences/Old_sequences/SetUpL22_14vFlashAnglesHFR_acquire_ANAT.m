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

P.numAcqs = 15;      % no. of Acquisitions in a Receive frame (this is a "superframe")
P.numFrames = 1;      % no. of Receive frames (real-time images are produced 1 per frame)
P.numAngle = 15;        % no. of flash angles for compounding

if (P.numAngle > 1)
    P.dtheta = (14*pi/180)/(P.numAngle-1);
    P.startAngle = -14*pi/180/2;
else
    P.dtheta = 0;
    P.startAngle=0;
end

P.prf=P.numAngle*500;  % pulse repetition frequency

P.maxCycle = 1; % Number superframes to collect
P.numCycle = 1;
P.duration = P.numAcqs/P.prf*P.maxCycle;
simulateMode = 0;   % set to acquire data using Vantage 128 hardware
visualize = 0; % 1 if want real time visualization

% -----------------------------------------------------------------------------

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L22-14v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 25;      % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';
% Media.MP(1,:) = [-45.5,0,30,1.0];
% Media.MP(2,:) = [-15.5,0,30,1.0];
% Media.MP(3,:) = [15.5,0,30,1.0];
% Media.MP(4,:) = [45.5,0,30,1.0];
% Media.NumPoints=4;

% Specify Resources.
Resource.Parameters.simulateMode = simulateMode;
Resource = setResources(Resource,P,visualize);
% Resource.InterBuffer(1).datatype = 'complex';
% Resource.InterBuffer(1).pagesPerFrame = P.numAcqs/P.numAngle;

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX = setTX(Resource,P,Trans);

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
                        'callMediaFunc', 1), 1,P.numFrames*P.numAcqs);  % movepoints EVERY acquisition to illustrate superframe concept
                                                                    % real-time images will look "jerky" but using the reconstructAll script,
                                                                    % playback process all acquisitions and shows smooth displacement

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
%     Receive(P.numAcqs*(i-1) + 1).callMediaFunc = 1;  % move points only once per frame
    for j = 1:P.numAcqs
        % -- Acquisitions for 'super' frame.
        rcvNum = P.numAcqs*(i-1) + j;
        Receive(rcvNum).Apod(:)=1;
        Receive(rcvNum).framenum = i;
        Receive(rcvNum).acqNum = j;
    end
end

Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:P.numAcqs);

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, P.numAcqs);
               
% % - Set specific ReconInfo attributes.
if P.numAngle>1
    
    for j = 1:P.numAcqs  % For each row in the column
        
        if isequal(rem(j,P.numAngle),1)
            ReconInfo(j).mode = 'replaceIQ'; % replace IQ data every P.numAngle acquisitions
        end
        
        if isequal(rem(j,P.numAngle),0) % cycle through angle acquisitions in each page
            ReconInfo(j).txnum=P.numAngle;
        else
            ReconInfo(j).txnum=rem(j,P.numAngle);
        end
        
        ReconInfo(j).rcvnum = j;
        ReconInfo(j).pagenum = ceil(j/P.numAngle);
        
        if isequal(rem(j,P.numAngle),0) % Acummulate angle data every na acquisitions
            ReconInfo(j).mode = 'accumIQ_replaceIntensity';
        end
        
    end
        
else
    ReconInfo(1).mode = 'replaceIntensity';
end

% Save Data Process
Process(1).classname = 'External';
Process(1).method = 'saveRFperFrame'; % calls the 'saveRFperFrame' function
Process(1).Parameters = {'srcbuffer','receive',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',-1,... % process the most recent frame.
    'dstbuffer','none'};

Process(2).classname = 'External';
Process(2).method = 'saveIQData';
Process(2).Parameters = {'srcbuffer','inter',... % name of buffer to process.
                    'srcbufnum',1,...
                    'srcframenum',1, ...
                    'dstbuffer','none'};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between acquisitions
SeqControl(2).argument = 1/P.prf*1e6;  %  usecs
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'sync';
SeqControl(4).argument = 1e6;   % 1s
nsc = 5; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numAcqs
        Event(n).info = 'Acquire RF';
        
        if isequal(rem(j,P.numAngle),0)
            Event(n).tx = P.numAngle;
        else
            Event(n).tx=rem(j,P.numAngle);
        end
        
        Event(n).rcv = P.numAcqs*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [2,nsc];
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        nsc = nsc + 1;

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 4;
    n = n+1;

    Event(n).info = 'Save RF data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n=n+1;

    % Do reconstruction and processing for 1st sub frame
    Event(n).info = 'Recon and save IQ data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 2;
    Event(n).seqControl = 3;
    n = n+1; 
end


EF(1).Function = text2cell('%saveRFperFrame');
EF(2).Function = text2cell('%saveIQData');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = P.numAcqs;

% Save all the structures to a .mat file.
% and invoke VSX automatically
filename = 'L22-14vFlashHFR_saveRFperFrame_acquire'; % used to launch VSX automatically
save(['MatFiles/',filename]);
VSX
commandwindow  % just makes the Command window active to show printout


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

%saveRFperFrame - save RF
saveRFperFrame(RData)

%Read in the misc variables struct 
endSample = evalin('base','Receive(end).endSample'); endSample

frameNum = 1;
foldername = ['Data/AnatData_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')];
mkdir(foldername)
RFfilename=[foldername '/RFdata_' foldername(15:end)];

P = evalin('base','P');
P.RFfilename=RFfilename;
P.foldername=foldername
assignin('base','P',P);

fname = ['frame_',num2str(frameNum,'%03.0f')];
RData=RData(1:endSample,:);
eval([fname,' = RData;']);

tic
save(RFfilename,fname,'-v6');

fprintf('saving time for RF frame %g = %g s \n',frameNum, toc);
return
%saveRFperFrame

%saveIQData
saveIQData(IQData)

figure(2); imagesc(abs(IQData))

frameNum = 1;
P = evalin('base','P');
foldername=P.foldername;
IQfilename=[foldername '/IQdata_' foldername(15:end)];
P.IQfilename=IQfilename;
assignin('base','P',P);

fname = ['frame_',num2str(frameNum,'%03.0f')];
eval([fname,' = IQData;']);

tic
save(IQfilename,fname,'-v6');
P = evalin('base','P');
datafilename=[foldername '/data_' foldername(15:end) '_P'];
save(datafilename,'P','-v6');

fprintf('saving time for IQ frame %g = %g s \n',frameNum, toc);
return

%saveIQData