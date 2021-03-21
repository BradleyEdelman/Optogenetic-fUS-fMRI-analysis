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
P.endDepth = 64;   % This should preferrably be a multiple of 128 samples.

P.numAngle = 7;        % no. of flash angles for compounding
P.numAcqs = P.numAngle*200;      % no. of Acquisitions in a Receive frame (this is a "superframe")
P.numFrames = 1;      % no. of Receive frames (real-time images are produced 1 per frame)

if (P.numAngle > 1)
    P.dtheta = 2*pi/180;
    P.startAngle = -6*pi/180;
else
    P.dtheta = 0;
    P.startAngle=0;
end

P.prf=P.numAngle*500;  % pulse repetition frequency

P.maxCycle = 200; % Number superframes to collect
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
% PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).PDelta = [1, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Resources.
Resource.Parameters.simulateMode = simulateMode;
Resource = setResources(Resource,P,visualize);
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).pagesPerFrame = P.numAcqs/P.numAngle;

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [15, 0.67, 6, 1];

% Specify TX structure array.
TX = setTX(Resource,P,Trans);

% Specify TGC Waveform structure.
TGC.CntrlPts = [330 560 780 1010 1023 1023 1023 1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.

BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];
     
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'InputFilter', BPF1, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1,P.numFrames*P.numAcqs);  % movepoints EVERY acquisition to illustrate superframe concept
                                                                    % real-time images will look "jerky" but using the reconstructAll script,
                                                                    % playback process all acquisitions and shows smooth displacement

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
%     Receive(P.numAcqs*(i-1) + 1).callMediaFunc = 1;  % move points only once per frame
    for j = 1:P.numAcqs
        % -- Acquisitions for 'super' frame.
        rcvNum = P.numAcqs*(i-1) + j;
%         Receive(rcvNum).Apod(:)=1;
        Receive(rcvNum).framenum = i;
        Receive(rcvNum).acqNum = j;
    end
end

Recon = struct('senscutoff', 0.75, ...
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
                   'pagenum', 1, ...
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
Process(2).method = 'filterIQ'; % calls the 'saveRFperFrame' function
Process(2).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',1,...
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

%     Event(n).info = 'Sync';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 0;
%     Event(n).seqControl = 4;
%     n = n+1;

    % Do reconstruction and processing for 1st sub frame
    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 2;
    Event(n).seqControl = [3,4];
    n = n+1; 
    
end


% --- If this last event is not included, the sequence stops after one pass, and enters "freeze" state
%     Pressing the freeze button runs the "one-shot" sequence one more time
%     For live acquisition in mode 0, simply comment out the 'if/end' statements and manually freeze and exit when the data looks good.
if simulateMode==2 || simulateMode==0 %  In live acquisiton or playback mode, run continuously, but run only once for all frames in simulation
    
    Event(n).info = 'Jump back to first event';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    
    if P.numCycle <= P.maxCycle
        Event(n).seqControl = 1;
    else
        Event(n).seqControl = 0;
    end
    P.numCycle = P.numCycle + 1;
end

EF(1).Function = text2cell('%saveRFperFrame');
EF(2).Function = text2cell('%filterIQ');

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

%saveRFperFrame - save RF
saveRFperFrame(RData)

persistent RFfilename frameNum
toc
%Read in the misc variables struct 
endSample = evalin('base','Receive(end).endSample');
if isempty(frameNum)
    frameNum = 1;
    RFfilename = ['D:/Data/RFdata_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')];
    P = evalin('base','P');
    P.RFfilename=RFfilename;
    assignin('base','P',P);
else
    frameNum = frameNum + 1;
end

fname = ['frame_',num2str(frameNum,'%03.0f')];
RData=RData(1:endSample,:);
eval([fname,' = RData;']);

tic
framelimit=300;
if isequal(frameNum,1)
    save(RFfilename,fname,'-v6');
elseif frameNum<framelimit
    save(RFfilename,fname,'-v6','-append');
elseif frameNum==framelimit
    save(RFfilename,fname,'-v6','-append');
    P = evalin('base','P');
    save([RFfilename '_P'],'P','-v6');
end
TOC=toc;
TIC=tic;
fprintf('saving time for frame %g = %g s \n',frameNum, toc);
return
%saveRFperFrame

%filterIQ
filterIQ(IQData)

persistent frameNum2 IQSAVE times TIC2
TOC=toc;
IQR=squeeze(IQData);

% figure(2); subplot(1,2,1);
% imagesc(abs(IQR(:,:,1))); 
% colormap gray
% title('Ultrasound image');
Tot=160;
if isempty(frameNum2)
    frameNum2 = 1;
    IQSAVE=zeros(size(IQR,1),size(IQR,2),size(IQR,3),Tot);
    times=zeros(1,Tot);
else
    frameNum2 = frameNum2 +1;
end

process=0; % 1 for filter, 2 for svd
if isequal(process,1)
    % computation of the power Doppler image
    [B,A]=butter(4,75/250,'high');    %coefficients for the high pass filter

    % sustraction of the first image
    % the signal stat at 0 and minimises filter oscilatons
    for i=1:size(IQR,1)
        for j=1:size(IQR,2)
            IQR(i,j,:)=IQR(i,j,:)-IQR(i,j,1); 
        end
    end

    sb=filter(B,A,IQR,[],3);    % blood signal (filtering in the time dimension)
    sb=sb(:,:,4:end);           % the first 4 temporal samples are eliminates (filter oscilations)
    PDI=mean(abs(sb).^2,3);     % computing the intensity of the blood signal the 

elseif isequal(process,2)
    
    [nz,nx,nt]=size(IQR);

    m=zeros(nx*nz,nt);
    i=0;
    for ix=1:nx
        for iz=1:nz
            i=i+1;
            m(i,:)=squeeze(IQR(iz,ix,:)-mean(IQR(iz,ix,:)));
        end
    end

    [U,S,V]=svd(m,'econ');
    
    nfilt=50;
    Sf=S;
    Sf(1:nfilt,1:nfilt)=0;
    mf=U*Sf*V';

    i=0;
    IQR1=zeros(nz,nx,nt);
    for ix=1:nx
        for iz=1:nz
            i=i+1;
            IQR1(iz,ix,:)=squeeze(mf(i,:));
        end
    end

    [B,A]=butter(4,0.3,'high');
    IQR2=filter(B,A,IQR1,[],3);
    IQR2=IQR2(:,:,4:end);           % the first 4 temporal samples are eliminates (filter oscilations)
    PDI=mean(abs(IQR2).^2,3);     % computing the intensity of the blood signal the 

end

% display the power doppler image
% subplot(1,2,2);
% figure;
% imagesc(10*log10(PDI./max(PDI(:)))); 
% % caxis([-35 0]);
% caxis auto
% colormap gray;
% title('PDI image')
if frameNum2<=Tot
    IQSAVE(:,:,:,frameNum2)=IQR;
end
TIC2=tic;
times(frameNum2)=TOC;
if frameNum2==Tot
    assignin('base','IQSAVE',IQSAVE);
    assignin('base','times',times);
end

if frameNum2==1
    fprintf('frame #%.f, Duration: %.2f, Cumulative: %.2f\n',frameNum2,times(1),times(frameNum2)-times(1)+mean(diff(times(1:frameNum2))));
else
    fprintf('frame #%.f, Duration: %.2f, Cumulative: %.2f\n',frameNum2,times(frameNum2)-times(frameNum2-1),times(frameNum2)-times(1)+mean(diff(times(1:frameNum2))));
end
%filterIQ






