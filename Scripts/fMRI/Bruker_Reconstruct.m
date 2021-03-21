path = '/home/bradley/Documents/fMRI/';
sessions = {'ofMRI_CW_Param_testing/20191010_093923_BEd_ofMRI_20191010_ofMRI_1_1_1';
    'ofMRI_CW_Param_testing/20191010_120101_BEd_ofMRI_20191010_ofMRI_2_1_2';
    'ofMRI_CW_Param_testing/20191010_135945_BEd_ofMRI_20191010_ofMRI_3_1_3';
    'ofMRI_CW_Param_testing/20191010_102736_BEd_ofMRI_20191010_ofMRI_5_1_5'};


% Determine file IDs to reconstruct
for  i = 1:size(sessions,1)
    
    raw_path = [path sessions{i} '/'];
    
    files = dir(raw_path);
    filenames = {files.name};
    folders = find(~isnan(str2double(filenames)) == 1);
    for j = size(folders,2):-1:1
        
        testpath = [raw_path filenames{folders(j)}];
        [DIM FovCm SliceThickMm VOX TR NEX] = BrukerInfo(testpath, 1);

        prefix = strcat(raw_path, filenames{folders(j)},'/pdata/', int2str(1), '/');
        loadname = sprintf('%s%s',prefix,'2dseq'); fid=fopen(loadname,'r'); imgtmp=fread(fid,'int16');
        if ~isequal(size(imgtmp,1),DIM(1)*DIM(2)*DIM(3)*NEX)
            folders(j) = [];
        end
    end    

    for k = 1:size(folders,2)
        Bruker2Nifti(raw_path,int2str(str2double(filenames(folders(k)))),1)
    end
    
end