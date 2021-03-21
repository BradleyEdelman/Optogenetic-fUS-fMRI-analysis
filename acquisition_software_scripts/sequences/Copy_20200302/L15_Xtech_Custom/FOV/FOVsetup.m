function Param = FOVsetup

% [filename, pathname] = uigetfile('D:\Data\*.mat', 'Select Anatomical Image for FOV Optimization');
% filenameParam = strcat(filename(1:end-4),'_P.mat');

% if isequal(filename,0)
    
%     % Specify PData structure array.
%     Param.delta = [.5, 0, 0.25];
%     Param.size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
%     Param.size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
%     Param.size(3) = 1;      % single image page
%     Param.origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
%     
% else
    
    % Check that file is an anatomical image dataset
    data = load(strcat(pathname,filename));
    frames = fieldnames(data);
    anat = data.(frames{end});
    if size(anat,1) < 100 || size(anat,2) < 100 || size(anat,3) < 50
        error('No Anatomy PDI Detected');
    else
        FOV({'pathname',pathname,'filename',filename,'filenameP',filenameParam})
    end

end






    
    