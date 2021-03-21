function Resource = setResources(Resource,P,visualize)

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;

% Specify Resources
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 18*P.endDepth*P.numAcqs;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;       % number of 'super frames'
Resource.InterBuffer(1).numFrames = 1;  % only one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = P.numFrames;

% Visualization window 
if isequal(visualize,1)
    Resource.DisplayWindow(1).Title = 'L22-14vFlashAnglesHFR';
    Resource.DisplayWindow(1).pdelta = 0.35;
    ScrnSize = get(0,'ScreenSize');
    DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
    DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
    Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                          DwWidth, DwHeight];
    Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
    Resource.DisplayWindow(1).Type = 'Verasonics';
    Resource.DisplayWindow(1).numFrames = 20;
    Resource.DisplayWindow(1).AxesUnits = 'mm';
    Resource.DisplayWindow(1).Colormap = gray(256);
end