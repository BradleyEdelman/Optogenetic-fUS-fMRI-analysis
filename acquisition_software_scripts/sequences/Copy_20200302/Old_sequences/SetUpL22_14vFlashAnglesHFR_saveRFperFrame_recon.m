% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name RunSetUpL11_5vFlashHFR_reconAll.m: "High Frame Rate" post-acquisition processing script (script 2 of 2)
%
% Description:
% To support higher acquisition frame rates with reduced DMA overhead, the first script (SetupL11_5vFlashHFR_acquire) acquires
%   numAcqs acquisitions into each RcvBuffer "super" frame, performing a transferToHost after each numAcqs set. After quitting VSX,
%   the RF data is saved in a matfile called 'RFdataHFR'. The real-time image reconstruction only processes the first subframe
%   in the super frame to minimize processing time and yet provide real-time feedback during data collection.
%   The data is saved in a matfile to permit the "clear all" at the beginning of every script. Note that in this HFR script,
%   each acquisition will lead to an image frame, and each DMA data transfer to the host contains one receive "super" frame.
%
% This second script has been written to reconstruct and displays all acquisitions. To be compatible with the acquired data,
%   it must use the same RcvBuffer size and is intended to be run in simulateMode=2 using the previously collected data. It works
%   in real-time as well, but the processing cannot keep up with acquisition, and some superframes are skipped but the script
%   displays all sub-frames within a super frame.
%
% To illustrate the HFR effect in simulation, the 'movePoints' function is invoked (in SetupL11_5vFlashHFR_acquire)
%   for every acquisition, resulting in a very discontinous apparent motion of the scatterers in "real-time" when processing
%   only one acquisition per superframe. However, the motion is smooth on playback when processing each T/R acquisition.
%
% last modification 12/13/2015

clear all

simulateMode = 2; % use mode=2 for post processing data collected by 'SetupL11_5vFlashHFR' and stored in 'RFdataHFR.mat'.
                  % use mode=0 for acquisiton, and processing of every acquisition until superframes are skipped

if simulateMode ==2
    
    folder='Data/RFdata_23-July-2018_15-54-59';
    numFile=1;
    file=[folder '/RFdataFlashAnglesHFR_' num2str(numFile)];
    load (file)  % this file includes RcvData and the numAcqs and numFrames parameters for properly setting Resources
    P.datafolder=folder;
    P.numFrames=size(RcvData{1,1},3);
    
else
    P.startDepth = 5;   % Acquisition depth in wavelengths
    P.endDepth = 64;   % This should preferrably be a multiple of 128 samples.
    P.numAcqs = 100;  % default no. of acquisitions in a frame (this is a "superframe") for testing
    P.numFrames = 10;  % default no. of frames
end

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
PData(1).PDelta = [.5, 0, 0.25];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Resources.
Resource.Parameters.simulateMode = simulateMode;
Resource = setResources(Resource,P,0);
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).pagesPerFrame = P.numAcqs/P.numAngle;

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [15,.67,6,1];

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
                        'callMediaFunc', 0), 1,P.numFrames*P.numAcqs);
                    
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numAcqs
        % -- Acquisitions for 'super' frame.
        rcvNum = P.numAcqs*(i-1) + j;
        Receive(rcvNum).Apod(:)=1;
        Receive(rcvNum).framenum = i;
        Receive(rcvNum).acqNum = j;
    end
end


% Specify Recon structure arrays: only need to define Recon for one Receive frame (superframe)
% - We need one Recon structure for each image to be produced, and here define all Recons needed for one superframe.
% In this Flash acquisition script, each acquisition is reconstructed and displayed as an image, and therefore
% each Recon will need just one ReconInfo structure, since every acquisition will be used to create a unique image.
% So, define ReconInfos in the same loop as Recons
% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
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



%  for acqNum = 1:P.numAcqs                   % Use this loop if independent Process structures are desired for explicit framenum
%      Process(acqNum) = Process(1);        % designation. The specific process structures must be referenced in the event sequence.
%      Process(acqNum).Parameters{4} = acqNum;
%  end

% Specify an external processing event.
Process(1).classname = 'External';
Process(1).method = 'saveIQData';
Process(1).Parameters = {'srcbuffer','inter',... % name of buffer to process.
                    'srcbufnum',1,...
                    'srcframenum',1, ...
                    'dstbuffer','none'};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 1/P.prf*1e6;  % usecs
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'sync';
SeqControl(4).argument = 1.1e6;   % 1s
nsc = 5; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:P.numFrames
    for j = 1:P.numAcqs  % Acquire all acquisitions defined in one super frame
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
        SeqControl(nsc).command = 'transferToHost'; % defines the receive frame
        nsc = nsc + 1;
        
    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 4;
    n = n+1;
    
    % Reconstruct
    Event(n).info = 'Recon and save IQ data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 3;
    
    
end

EF(1).Function = text2cell('%saveIQData');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = P.numAcqs;

% Save all the structures to a .mat file.
filename = 'L22-14vFlashHFR_saveRFperFrame_recon'; % used to launch VSX automatically
save(['MatFiles/',filename]);
disp([ mfilename ': NOTE -- Running VSX automatically!'])
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

%saveIQData
saveIQData(IQData)
    
persistent IQfilename frameNum

if isempty(frameNum)
    frameNum = 1;
    P = evalin('base','P');
    RFfilename=P.RFfilename;
    IQfilename=[P.datafolder '/IQ' RFfilename(8:end)]
    P.IQfilename=IQfilename;
    assignin('base','P',P);
else
    frameNum = frameNum + 1;
end

fname = ['frame_',num2str(frameNum,'%03.0f')];
eval([fname,' = IQData;']);

tic
framelimit=5;
if isequal(frameNum,1)
    save(IQfilename,fname,'-v6');
elseif frameNum<framelimit
    save(IQfilename,fname,'-v6','-append');
elseif frameNum==framelimit
    save(IQfilename,fname,'-v6','-append');
    P = evalin('base','P');
    save([IQfilename '_P'],'P','-v6');
end
fprintf('saving time for frame %g = %g s \n',frameNum, toc);
return

%saveIQData