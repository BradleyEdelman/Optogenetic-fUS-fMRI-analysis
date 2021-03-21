function [DIM FovCm SliceThickMm VOX TR NEX] = BrukerInfo(path, recofolder)

% get commonly used parameters of protocol from Brucker Raw data folder
% 2009/03/27, by Iris

% path
strcat(path,'/pdata/', int2str(recofolder))
reco = strcat(path,'/pdata/', int2str(recofolder), '/','reco')
file = textread(reco,'%s','delimiter','\n','whitespace','');

for i=1:length(file)
    tmp = cell2mat(file(i));

    % DIM: 
    if length(tmp)>=12 & findstr(tmp, '##$RECO_size');
            xy=str2num(cell2mat(file(i+1)));
    end   

end


fn = strcat(path,'/method');
file = textread(fn,'%s','delimiter','\n','whitespace','');
TR=0; 
for i=1:length(file)
    tmp = cell2mat(file(i));

%% FovCm
     if length(tmp)>=12 & findstr(tmp, '##$PVM_FovCm');
            fov_xy=str2num(cell2mat(file(i+1)));
     end 
     
     %% input matrix size
     if length(tmp)>=12 & findstr(tmp, '##$PVM_Matrix');
            xy_input=str2num(cell2mat(file(i+1)));
     end 
     
     %% TR
     
     if length(tmp)>=10 & findstr(tmp, '##$PVM_RepetitionTime');
           pos1=findstr(tmp, '=');            % pulse excitation TR, smaller than scan TR in 3D; equals to scan TR in 2D
           TR= str2num(tmp(pos1+1:end))/1000; % in terms of second
     end      
     
     if length(tmp)>=10 & findstr(tmp, '##$Scan_RepetitionTime');
           pos1=findstr(tmp, '=');
           TR= str2num(tmp(pos1+1:end))/1000; % in terms of second
     end 
     
    %% Slice Thickness   
    if length(tmp)>=17 & findstr(tmp, '##$PVM_SliceThick');
           pos1=findstr(tmp, '=');
           posdot=findstr(tmp(pos1+1:end),'.');

           if isempty(posdot)
               SliceThickMm= str2num(tmp(pos1+1:end));
           else
               SliceThickMm = str2num(tmp(pos1+1:posdot-1))+str2num(tmp(posdot+1:end))/10^size(tmp(posdot+1:end),2);
           end
    end
    
    %% Slice distance
     if length(tmp)>=22 & findstr(tmp, '##$PVM_SPackArrSliceDistance');
        tmp_slicedis = str2num(cell2mat(file(i+1)));
        if size(tmp_slicedis,2)==1
            vox_z=tmp_slicedis;
        else
            vox_z=sum(tmp_slicedis);
        end
     end
    
    %% No. of slices
    if length(tmp)>=22 & findstr(tmp, '##$PVM_SPackArrNSlices');
        tmp_slice = str2num(cell2mat(file(i+1)));
        if size(tmp_slice,2)==1
            DIM(3)=tmp_slice;
        else
            DIM(3)=sum(tmp_slice);
        end
        
%         if DIM(3)==1 && (size(tmp_dim,2)>2)
%             DIM(3) = tmp_dim(3);
%             tmp=DIM(1:2);
%             DIM(1:2)=tmp(2:-1:1);% if bruker 3D files' first 2 dimension is (256*128),then data in matlab is the opposite(128*256) 
%         end
    end
    
    %% No. of timepoints
    if length(tmp)>=19 & findstr(tmp, '##$PVM_NRepetitions');
           pos1=findstr(tmp, '=');
           NEX = str2num(tmp(pos1+1:end));
    else
        NEX=1;
    end
    
     %% No. of average
    if length(tmp)>=17 & findstr(tmp, '##$PVM_NAverages');
           pos1=findstr(tmp, '=');
           Aver = str2num(tmp(pos1+1:end));
    end
    
 %% Mode
    if length(tmp)>=18 & findstr(tmp, '##$PVM_SpatDimEnum');
           pos1=findstr(tmp, '=');
           Mode = tmp(pos1+1:end);
    end   
   
end


%%  Mode 2D or 3D
if Mode == '<2D>'
    DIM=[xy tmp_slice];%xy=RECO size
    FovCm=[fov_xy vox_z*tmp_slice/10];
    VOX = [fov_xy*10./xy vox_z];% unit in mm
    TR=TR*Aver;%TrueFISP
else  %3D
    DIM=[xy];
    FovCm=[fov_xy];
    VOX = [fov_xy*10./xy];% unit in mm
    TR=TR*xy_input(3)*Aver;%3D TrueFISP
end

%%

