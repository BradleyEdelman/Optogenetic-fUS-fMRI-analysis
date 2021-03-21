function [handles, hObject] = FOVinit(handles, hObject)

if isequal(handles.data.pathname,0) || isequal(handles.data.filename,0)
    fprintf('\nNo file selected\n')
else
    % Load and store image
    fileAnat = strcat(handles.data.pathname, handles.data.filename);
    % Check data file quickly
    data = load(fileAnat);
    frames = fieldnames(data);
    anat = data.(frames{end});
    if size(anat,1) < 100 || size(anat,2) < 100 || size(anat,3) < 50
		fprintf('No Anatomy PDI Detected');
    else
        [pdi,times,block,PSD,velocity] = Clutterfilt(fileAnat,'filtcutoff',0.3,'velocity',0);
        pdi=10*log10(pdi./max(pdi(:)));
        axes(handles.axes1); axis off
        handles.data.image = imagesc(pdi); colormap gray; caxis([-25 0])
        handles.data.pdi = pdi;

        % Load and set parameters
        fileParam = strcat(handles.data.pathname, handles.data.filenameP);
        load(fileParam)

        % Depth specification
        handles.data.depthStart = P.startDepth;
        set(handles.depthstart,'string',num2str(P.startDepth));
        handles.data.depthEnd = P.endDepth;
        set(handles.depthend,'string',num2str(P.endDepth));

        % Spatial sampling delta
        handles.data.xdelta = .5; % P.xdelta;
        set(handles.xdelta,'string',handles.data.xdelta);
        handles.data.ydelta = 0; % P.ydelta;
        set(handles.ydelta,'string',handles.data.ydelta);
        handles.data.zdelta = .25; % P.zdelta;
        set(handles.zdelta,'string',handles.data.zdelta);

        handles.data.xdeltatmp = handles.data.xdelta;
        handles.data.ydeltatmp = handles.data.ydelta;
        handles.data.zdeltatmp = handles.data.zdelta;

        % Initial FOV size
        [h, w] = size(pdi);
        set(handles.xFOV,'string',num2str(w));
        set(handles.zFOV,'string',num2str(h));

        % Set parameters of sliders
        set(handles.sliderZtop,'min',h/2,'max',h,'sliderstep',[1/(h-h/2) 10/(h-h/2)],'value',h);
        handles.data.lineZtop = line([1 w],[1 1],'color','r');

        set(handles.sliderZbottom,'min',1,'max',h/2,'sliderstep',[1/(h/2-1) 10/(h/2-1)],'value',1);
        handles.data.lineZbottom = line([1 w],[h h],'color','r');

        set(handles.sliderXleft,'min',1,'max',w/2,'sliderstep',[1/(w/2-1) 10/(w/2-1)],'value',1);
        handles.data.lineXleft = line([1 1],[1 h],'color','r');

        set(handles.sliderXright,'min',w/2,'max',w,'sliderstep',[1/(w-w/2) 10/(w-w/2)],'value',w);
        handles.data.lineXright = line([w w],[1 h],'color','r');

    end
end


