if isunix
    cd('/media/bradley/Seagate Backup Plus Drive/fMRI/')
    addpath(genpath('/media/bradley/Seagate Backup Plus Drive/Preprocessing/'))
    slash = '/';
elseif ispc
    cd('D:\fMRI\');
    addpath(genpath('D:\Preprocessing'));
    slash = '\';
end

data_fold = ['ofMRI_pre_CW' slash];
% data_fold = [ 'ofMRI_post_CW' slash];

sessions = {'20200211_090028_BEd_preCW_4419410_N_1_1_64';
    '20200210_161222_BEd_preCW_4419409_N_1_1_63';
    };

for  j = 1:size(sessions,1)
    
    raw  = dir([pwd slash data_fold sessions{j} slash]);
    filenames = {raw.name};
    folders = find(~isnan(str2double(filenames)) == 1);

    for i = 1:size(folders,2)

        inputfolder = [pwd slash data_fold sessions{j} slash filenames{folders(i)}];
        Read_Bruker(inputfolder)

    end

end










