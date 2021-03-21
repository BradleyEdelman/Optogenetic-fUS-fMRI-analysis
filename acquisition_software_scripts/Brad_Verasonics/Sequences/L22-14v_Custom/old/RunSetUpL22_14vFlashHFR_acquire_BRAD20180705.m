% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: RunSetUpL11_5vFlashHFR_acquire.m - "High Frame Rate" acquisition for ultrafast imaging (script 1 of 2).
%
% Description:
% NOTE: This is a unique example that runs VSX automatically, and after closing the Control GUI, runs the post-processing
%       script (script 2 of 2) automatically. Two scripts are part of this example!!!
%
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
% A second script for reconstructing and displaying all aquisitions has been developed to be run in simulateMode=2 after this one,
%   using the previously collected data. That script, SetUpL11_5vFlashHFR_reconAll.m, reconstructs and displays all
%   acquisitions by defining ReconInfo structures for each acquisition, and defining a Recon event for each display frame.
%   Acquiring data using that script results in skipped superframes when processing cannot keep up with acquisition.
%
% To illustrate the HFR acquisition and display in simulation, the 'movePoints' function is invoked for every acquisition,
%   resulting in a very discontinous apparent motion of the scatterers in "real-time".
%   However, the motion is smooth on playback when post-processing each T/R acquisition with the "reconAll" script.
%
% For convenience, this script is currently set to launch VSX automatically. Upon quitting VSX (by closing the control GUI), the "reconAll"
%   script is automatically invoked to process and display all of the frames just collected.
%   Note that in simulateMode=1, this script runs once, filling all numFrames frames, then freezing. Simply close the GUI, and
%   the reconAll script launches automatically to demonstrate the post-processing of all acquisitions.
%
% Last modification 12/13/2015

clear all

% --- Frequently Modified Parameters ------------------------------------------
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 128;   % This should preferrably be a multiple of 128 samples.

P.numAcqs = 1000;      % no. of Acquisitions in a Receive frame (this is a "superframe")
P.numFrames = 300;      % no. of Receive frames (real-time images are produced 1 per frame)
P.numAngle = 5;        % no. of flash angles for compounding

if (P.numAngle > 1)
    P.dtheta = (12*pi/180)/(P.numAngle-1);
    P.startAngle = -12*pi/180/2;
else
    P.dtheta = 0;
    P.startAngle=0;
end

P.prf=P.numAngle*500;  % pulse repetition frequency

% P.maxCycles = 2; % Number superframes to collect
P.duration = P.numAcqs/P.prf*P.numFrames;
simulateMode = 0;   % set to acquire data using Vantage 128 hardware
visualize = 0; % 1 if want real time visualization

RFdataFilename = 'RFdataHFR';
if simulateMode==2      % playback using this script still results in processing only one acquisiton per super frame
    load (RFdataFilename)
end

% -----------------------------------------------------------------------------

filename = mfilename; % used to launch VSX automatically

% Specify Trans structure array.
Trans.name = 'L22-14v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % L22-14v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 25;      % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Resources.
Resource.Parameters.simulateMode = simulateMode;
Resource = setResources(Resource,P,visualize);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX = setTX(Resource,P,Trans);

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,298,395,489,618,727,921,1023];
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
    % Receive(P.numAcqs*(i-1) + 1).callMediaFunc = 1;  % move points only once per frame
    for j = 1:P.numAcqs
        % -- Acquisitions for 'super' frame.
        rcvNum = P.numAcqs*(i-1) + j;
        Receive(rcvNum).Apod(:)=1;
        Receive(rcvNum).framenum = i;
        Receive(rcvNum).acqNum = j;
    end
end

% Specify Recon structure arrays.
if isequal(visualize,1)
    Recon = struct('senscutoff', 0.6, ...
                   'pdatanum', 1, ...
                   'rcvBufFrame', -1, ...     % use most recently transferred frame
                   'IntBufDest', [1,1], ...
                   'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
                   'RINums', 1);

    Define ReconInfo structures.  (just one ReconInfo to process one acquisition per superframe)
    ReconInfo(1) = struct('mode', 'replaceIntensity', ...
                       'txnum', 1, ...
                       'rcvnum', 1, ...  % use the first acquisition of each frame
                       'regionnum', 1);
end

% Specify Process structure array.
if isequal(visualize,1)
    pers = 20;
    Process(1).classname = 'Image';
    Process(1).method = 'imageDisplay';
    Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                             'framenum',-1,...   % (-1 => lastFrame)
                             'pdatanum',1,...    % number of PData structure to use
                             'pgain',1.0,...            % pgain is image processing gain
                             'reject',2,...      % reject level
                             'persistMethod','simple',...
                             'persistLevel',pers,...
                             'interpMethod','4pt',...
                             'grainRemoval','none',...
                             'processMethod','none',...
                             'averageMethod','none',...
                             'compressMethod','power',...
                             'compressFactor',40,...
                             'mappingMode','full',...
                             'display',1,...      % display image after processing
                             'displayWindow',1};
end

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between acquisitions (& superframes)
SeqControl(2).argument = 200;  % 200 usecs
SeqControl(3).command = 'triggerOut';
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

n = 1; % n is count of Events
% Cycles=1;
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
        Event(n).seqControl = [2,3];
        n = n+1;
        
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [2,3,nsc];
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        nsc = nsc + 1;
    
    % Do reconstruction and processing for 1st sub frame
    if isequal(visualize,1)
        Event(n).info = 'Reconstruct';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 1;
        Event(n).process = 1;
        Event(n).seqControl = 4;
        n = n+1;
    end
%     Cycles = Cycles + 1;
end

% --- If this last event is not included, the sequence stops after one pass, and enters "freeze" state
%     Pressing the freeze button runs the "one-shot" sequence one more time
%     For live acquisition in mode 0, simply comment out the 'if/end' statements and manually freeze and exit when the data looks good.
% if simulateMode==2 || simulateMode==0 %  In live acquisiton or playback mode, run continuously, but run only once for all frames in simulation
%     Event(n).info = 'Jump back to first event';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 0;
%     Event(n).seqControl = 1;
% end

% User specified UI Control Elements
% - Sensitivity Cutoff
if isequal(visualize,1)
    UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                      'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                      'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
    UI(1).Callback = text2cell('%SensCutoffCallback');

%     - Range Change
    MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
    AxesUnit = 'wls';
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
            AxesUnit = 'mm';
            MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
        end
    end
    UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                     'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
    UI(2).Callback = text2cell('%RangeChangeCallback');
end

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = P.numAcqs;

% Save all the structures to a .mat file.
% and invoke VSX automatically
save(['MatFiles/',filename]);
disp([ mfilename ': NOTE -- Running VSX automatically!']), disp(' ')
VSX
commandwindow  % just makes the Command window active to show printout

 % after quitting VSX, save the RFData collected at high frame rate
% if simulateMode ~= 2
%     disp ('Info:  Saving the RF Data buffer -- please wait!'), disp(' ')
tic
    save (['MatFiles/',RFdataFilename], 'RcvData', 'P', '-v7.3');toc
% end

% then run full playback automatically after quitting the first script
% RunSetUpL22_14vFlashHFR_reconAll_BRAD20180705    % be sure that simulateMode=2 in that script!!!

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
