function varargout = FOV(varargin)
% FOV MATLAB code for FOV.fig
%      FOV, by itself, creates a new FOV or raises the existing
%      singleton*.
%
%      H = FOV returns the handle to a new FOV or the handle to
%      the existing singleton*.
%
%      FOV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOV.M with the given input arguments.
%
%      FOV('Property','Value',...) creates a new FOV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FOV_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FOV_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FOV

% Last Modified by GUIDE v2.5 18-Mar-2019 18:32:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FOV_OpeningFcn, ...
                   'gui_OutputFcn',  @FOV_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FOV is made visible.
function FOV_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FOV (see VARARGIN)

% Choose default command line output for FOV
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FOV wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FOV_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on slider movement.
function sliderXleft_Callback(hObject, eventdata, handles)
% hObject    handle to sliderXleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
[handles, hObject] = FOVupdate(handles, hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sliderXleft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderXleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLIDERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function sliderXright_Callback(hObject, eventdata, handles)
% hObject    handle to sliderXright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
[handles, hObject] = FOVupdate(handles, hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sliderXright_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderXright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderZtop_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZtop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
[handles, hObject] = FOVupdate(handles, hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sliderZtop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZtop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderZbottom_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
[handles, hObject] = FOVupdate(handles, hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sliderZbottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function depthstart_Callback(hObject, eventdata, handles)
% hObject    handle to depthstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depthstart as text
%        str2double(get(hObject,'String')) returns contents of depthstart as a double
depthstart = str2double(get(handles.depthstart,'string'));
if depthstart < handles.data.depthStart
    set(handles.hObject,'string',num2str(handles.data.depthStart));
else
    zdelta = str2double(get(handles.zdelta,'string'));
    depthstart = round(depthstart/zdelta)*zdelta;
    ztop = (handles.data.depthEnd - depthstart) / zdelta;
    set(handles.sliderZtop,'value',ztop)
end
[handles, hObject] = FOVupdate(handles, hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function depthstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depthstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function depthend_Callback(hObject, eventdata, handles)
% hObject    handle to depthend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depthend as text
%        str2double(get(hObject,'String')) returns contents of depthend as a double
depthend = str2double(get(handles.depthend,'string'));
if depthend > handles.data.depthEnd
    set(handles.hObject,'string',num2str(handles.data.depthEnd));
else

    xdelta = round(str2double(get(handles.xdelta,'string'))*100)/100;
    xdsample = xdelta/handles.data.xdelta;
    
    zdelta = round(str2double(get(handles.zdelta,'string'))*100)/100;
    zdsample = zdelta/handles.data.zdelta;

    [h, w] = size(handles.data.pdi(1:zdsample:end,1:xdsample:end));
    zbottom = h + 1 - ((depthend - handles.data.depthStart) / zdelta);
    set(handles.sliderZbottom,'value',zbottom)
    
end
[handles, hObject] = FOVupdate(handles, hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function depthend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depthend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xFOV_Callback(hObject, eventdata, handles)
% hObject    handle to xFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xFOV as text
%        str2double(get(hObject,'String')) returns contents of xFOV as a double
xleft = get(handles.sliderXleft,'value');
xright = get(handles.sliderXright,'value');
set(handles.xFOV, 'string', xright - xleft + 1)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xFOV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function zFOV_Callback(hObject, eventdata, handles)
% hObject    handle to zFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zFOV as text
%        str2double(get(hObject,'String')) returns contents of zFOV as a double
zdelta = str2double(get(handles.zdelta,'string'));
depthend = str2double(get(handles.depthend,'string'));
depthstart = str2double(get(handles.depthstart,'string'));
set(handles.zFOV, 'string', (depthend - depthstart) / zdelta)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function zFOV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clipFOV.
function clipFOV_Callback(hObject, eventdata, handles)
% hObject    handle to clipFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of clipFOV
[handles, hObject] = FOVupdate(handles, hObject);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELTAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdelta_Callback(hObject, eventdata, handles)
% hObject    handle to xdelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xdelta as text
%        str2double(get(hObject,'String')) returns contents of xdelta as a double
[handles, hObject] = FOVupdate(handles, hObject);
handles.data.xdeltatmp = str2double(get(handles.xdelta,'string'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xdelta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xdelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ydelta_Callback(hObject, eventdata, handles)
% hObject    handle to ydelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ydelta as text
%        str2double(get(hObject,'String')) returns contents of ydelta as a double

% --- Executes during object creation, after setting all properties.
function ydelta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ydelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function zdelta_Callback(hObject, eventdata, handles)
% hObject    handle to zdelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zdelta as text
%        str2double(get(hObject,'String')) returns contents of zdelta as a double
[handles, hObject] = FOVupdate(handles, hObject);
handles.data.zdeltatmp = str2double(get(handles.zdelta,'string'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function zdelta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zdelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD/SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in loadPDI.
function loadPDI_Callback(hObject, eventdata, handles)
% hObject    handle to loadPDI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('D:\Data\*.mat', 'Select Anatomical Image for FOV Optimization');
filenameParam = strcat(filename(1:end-4),'_P.mat');
handles.data.pathname = pathname;
handles.data.filename = filename;
handles.data.filenameP = filenameParam;
[handles, hObject] = FOVinit(handles, hObject);
guidata(hObject, handles);


% --- Executes on button press in loadFOV.
function loadFOV_Callback(hObject, eventdata, handles)
% hObject    handle to loadFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('D:\Data\*FOV.mat', 'Select Anatomical Image for FOV Optimization');
load(strcat(pathname,filename));

if strcmp(pathname,P.localizerpathname)
    handles.data.pathname = P.localizerpathname;
elseif contains(pathname,P.localizerpathname)
    handles.data.pathname = pathname;
end
handles.data.filename = P.localizerfilename;
handles.data.filenameP = P.localizerfilenameP;

[handles, hObject] = FOVinit(handles, hObject);

set(handles.depthstart,'string',num2str(P.depthstart));
depthstart_Callback(hObject, eventdata, handles)

set(handles.depthend,'string',num2str(P.depthend));
depthend_Callback(hObject, eventdata, handles)

set(handles.sliderXleft,'value',P.sliderX(1)*(P.delta(1)/handles.data.xdeltatmp));
sliderXleft_Callback(hObject, eventdata, handles)

set(handles.sliderXright,'value',P.sliderX(2)*(P.delta(1)/handles.data.xdeltatmp));
sliderXright_Callback(hObject, eventdata, handles)

set(handles.zdelta,'string',num2str(P.delta(3)));
zdelta_Callback(hObject, eventdata, handles)
handles.data.zdeltatmp = str2double(get(handles.zdelta,'string'));

set(handles.ydelta,'string',num2str(P.delta(2)));
ydelta_Callback(hObject, eventdata, handles)

set(handles.xdelta,'string',num2str(P.delta(1)));
xdelta_Callback(hObject, eventdata, handles)
handles.data.xdeltatmp = str2double(get(handles.xdelta,'string'));

guidata(hObject, handles);


% --- Executes on button press in saveFOV.
function saveFOV_Callback(hObject, eventdata, handles)
% hObject    handle to saveFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Sampling info
P.delta = [str2double(get(handles.xdelta,'string')) ...
    str2double(get(handles.ydelta,'string'))...
    str2double(get(handles.zdelta,'string'))];

xdelta = round(str2double(get(handles.xdelta,'string'))*100)/100;
set(handles.xdelta,'string',num2str(xdelta))
xdsample = xdelta/handles.data.xdelta;
zdelta = round(str2double(get(handles.zdelta,'string'))*100)/100;
set(handles.zdelta,'string',num2str(zdelta));
zdsample = zdelta/handles.data.zdelta;
P.dsample = [xdsample 1 zdsample];


% Depth info
P.depthstart = str2double(get(handles.depthstart,'string'));
P.depthend = str2double(get(handles.depthend,'string'));
zbottomlim = [get(handles.sliderZbottom,'min') get(handles.sliderZbottom,'max')];
ztoplim = [get(handles.sliderZtop,'min') get(handles.sliderZtop,'max')];
zbottomval = round(get(handles.sliderZbottom,'value'));
ztopval = round(get(handles.sliderZtop,'value'));
P.sliderZ = [zbottomval ztopval];
P.sliderZbottomlim = zbottomlim;
P.sliderZtoplim = ztoplim;

% Width info
xleftlim = [get(handles.sliderXleft,'min') get(handles.sliderXleft,'max')];
xrightlim = [get(handles.sliderXright,'min') get(handles.sliderXright,'max')];
xleftval = round(get(handles.sliderXleft,'value'));
xrightval = round(get(handles.sliderXright,'value'));
P.sliderX = [xleftval xrightval];
P.sliderXleftlim = xleftlim;
P.sliderXrightlim = xrightlim;

% Origin info
P.origin = [-(xleftlim(2) - xleftval)*P.delta(1) 0 P.depthstart];
P.size(1) = ceil((P.depthend-P.depthstart)/P.delta(3));
P.size(2) = xrightval - xleftval + 1;
P.size(3) = 1;

% Image info
P.img.orig = handles.data.pdi;
imgtmp = handles.data.pdi;
[h, w] = size(imgtmp);
imgtmp(1:h-ztopval+1,:) = 0;
imgtmp(h-zbottomval+1:h,:) = 0;
imgtmp(:,1:xleftval) = 0;
imgtmp(:,xrightval:w) = 0;
P.img.clip = imgtmp;
imgtmp = handles.data.pdi;
[h, w] = size(imgtmp);
imgtmp(h-zbottomval+1:h,:) = [];
imgtmp(1:h-ztopval+1,:) = [];
imgtmp(:,xrightval:w) = [];
imgtmp(:,1:xleftval) = [];
P.img.cut = imgtmp;

P.localizerpathname = handles.data.pathname;
P.localizerfilename = handles.data.filename;
P.localizerfilenameP = handles.data.filenameP;
assignin('base','Param',P);

try
    Date = evalin('base','Date');
catch
    Date = datestr(now,'dd-mmmm-yyyy_HH-MM-SS');
end

handles.data.filenameFOV = strcat('IQData_',Date,'_FOV.mat');
save(strcat(handles.data.pathname,handles.data.filenameFOV),'P','-v7.3');
guidata(hObject, handles);

