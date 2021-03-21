function [handles, hObject] = FOVupdate(handles, hObject)

xdelta = round(str2double(get(handles.xdelta,'string'))*100)/100;
set(handles.xdelta,'string',num2str(xdelta))
xdsample = xdelta/handles.data.xdelta;

zdelta = round(str2double(get(handles.zdelta,'string'))*100)/100;
set(handles.zdelta,'string',num2str(zdelta));
zdsample = zdelta/handles.data.zdelta;

[h, w] = size(handles.data.pdi(1:zdsample:end,1:xdsample:end));
imgtmp = handles.data.pdi(1:zdsample:end,1:xdsample:end);

% Set parameters of sliders
ztop = round(get(handles.sliderZtop,'value')) / (zdelta/handles.data.zdeltatmp);
if ztop > h; ztop = h; end
set(handles.sliderZtop,'min',h/2,'max',h,'sliderstep',[1/(h-h/2) 10/(h-h/2)],'value',ztop);

zbottom = round(get(handles.sliderZbottom,'value')) / (zdelta/handles.data.zdeltatmp);
if zbottom < 1; zbottom = 1; end
set(handles.sliderZbottom,'min',1,'max',h/2,'sliderstep',[1/(h/2-1) 10/(h/2-1)],'value',zbottom);


depthstart = handles.data.depthEnd - (ztop * zdelta);
depthend = handles.data.depthStart + ((h-zbottom+1) * zdelta);
set(handles.depthstart,'string',num2str(depthstart));
set(handles.depthend,'string',num2str(depthend));
set(handles.zFOV, 'string', (depthend - depthstart) / zdelta)


xleft = round(get(handles.sliderXleft,'value')) / (xdelta/handles.data.xdeltatmp);
if xleft < 1; xleft = 1; end
set(handles.sliderXleft,'min',1,'max',w/2,'sliderstep',[1/(w/2-1) 10/(w/2-1)],'value',xleft);

xright = round(get(handles.sliderXright,'value')) / (xdelta/handles.data.xdeltatmp);
if xright > w; xright = w; end
set(handles.sliderXright,'min',w/2,'max',w,'sliderstep',[1/(w-w/2) 10/(w-w/2)],'value',xright);


set(handles.xFOV, 'string', xright - xleft + 1)


clipFOV = get(handles.clipFOV,'value');
if isequal(clipFOV,1)

    ztopimg = round(get(handles.sliderZtop,'value'));
    imgtmp(1:h-ztopimg+1,:) = 0;

    zbottomimg = round(get(handles.sliderZbottom,'value'));
    imgtmp(h-zbottomimg+1:h,:) = 0;

    xleftimg = round(get(handles.sliderXleft,'value'));
    imgtmp(:,1:xleftimg) = 0;

    xrightimg = round(get(handles.sliderXright,'value'));
    imgtmp(:,xrightimg:w) = 0;
    
end
handles.data.image = imagesc(imgtmp);

handles.data.lineZtop = line([1 w],[h-ztop+1 h-ztop+1],'color','r');
handles.data.lineZbottom = line([1 w],[h-zbottom+1 h-zbottom+1],'color','r');
handles.data.lineXleft = line([xleft xleft],[1 h],'color','r');
handles.data.lineXright = line([xright xright],[1 h],'color','r');







