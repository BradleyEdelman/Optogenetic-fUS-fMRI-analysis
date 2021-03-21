function Bruker2Nifti(prepath, ScanSession, recofolder)
%recofolder 
%Bruker2Nifti('E:\fMRI\Data\2012\20120625_3DbSSFP_protocal_test\Iris201207.eK2\29',2)

path=[prepath ScanSession]
%%Read format information
[DIM FovCm SliceThickMm VOX TR NEX] = BrukerInfo(path, recofolder)

%% Convert Bruker raw data to analyze format
prefix = strcat(prepath, ScanSession,'/pdata/', int2str(recofolder), '/');

% cd (prefix);
% filename='2dseq';

origin= [0 0 0];
datatype = 4;

loadname = sprintf('%s%s',prefix,'2dseq');
fid=fopen(loadname,'r');
imgtmp=fread(fid,'int16');

img1=reshape(imgtmp,DIM(1),DIM(2),DIM(3),NEX);
% img1=img1(:,:,:,1:380);


% for slice = 1:DIM(3)
%     for time=1:time_points
%         img1(:,:,slice,time) = fliplr(img1(:,:,slice,time));
% %         img1(:,:,slice,time) = rot90(img1(:,:,slice,time),2);
% %         img1(:,:,slice,time) = fliplr(img1(:,:,slice,time));
%         
%     end
% end

fclose(fid);
mkdir ([prepath 'Results/' ScanSession])
filename=[prepath 'Results/' ScanSession '/2dseq'];
nii = make_nii(img1,VOX,origin, datatype);
save_nii(nii, filename );
clear img1;
delete([filename '.mat']);

