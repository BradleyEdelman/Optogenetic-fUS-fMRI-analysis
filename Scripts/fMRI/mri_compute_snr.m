%%
clear all
% Main data directory
if isunix
    Dir = '/media/bradley/Seagate Backup Plus Drive/';
    slash = '/';
elseif ispc
    Dir = 'D:\';
    slash = '\';
end
cd(Dir)
addpath(genpath([Dir 'Preprocessing' slash]));
storage = [Dir 'Data_Processed' slash 'fMRI' slash]; if ~exist(storage,'dir'); mkdir(storage); end
% fmri_fold = [Dir 'fMRI' slash 'ofMRI_pre_CW' slash];
addpath(genpath(storage));

% Total "pre" ofMRI
base_fold = {'20191108_083343_BEd_preCW_4346075_N_D5_1_47';
    '20191108_143147_BEd_preCW_4364124_R1_D1_1_51';
    '20191108_162123_BEd_preCW_4364124_L1_D1_1_52';
    '20191111_151047_BEd_preCW_4364122_N_D2_1_53';
    '20191111_220250_BEd_preCW_4364123_R1_D2_2_1_57';
    '20191116_093138_BEd_preCW_4364121_N_D2_1_58';
    '20191116_133812_BEd_preCW_4364122_R1_D3_1_60';
    };

% Total "post" ofMRI
base_fold = {'20191122_095408_BEd_postCW_4346075_N_1_1';
    '20191122_123846_BEd_postCW_4364143_N_1_2';
    '20191122_143826_BEd_postCW_4364124_L1_1_3';
    '20191124_191348_BEd_postCW_4364122_N_1_4';
    '20191124_204714_BEd_postCW_4364123_R1_1_5';
    };

SNR = 1;



%% View all functional slices
for i = 1:size(scannum,1)
    
    for j = 1:size(scannum{i},2)
    
        basefile = [raw_data_fold{1} num2str(scannum{i}(j)) 'f\' num2str(scannum{i}(j)) '_EPI'];
        nii{i,j} = load_nii([basefile '.nii']);
        param{i,j} = load([basefile '_params.mat']);
        figure(1);
        for k = 1:size(nii{i,j}.img,3)
            subplot(5,ceil(size(nii{i,j}.img,3)/5),k)
            imagesc(nii{i,j}.img(:,:,k,1)')
            set(gca,'xticklabel',[],'yticklabel',[])
        end
        colormap gray
    
        idx(i,j) = input('\npick slice of interest\n');
    
        f = figure(2);
        imagesc(nii{i,j}.img(:,:,idx(i,j),1)'); axis tight; colormap gray
        [xsig{i,j}, ysig{i,j}] = getpts(f);
        [xnoise{i,j}, ynoise{i,j}] = getpts(f);

        t = linspace(0, 2*pi); r = 2;

        xsig{i,j} = r * cos(t) + xsig{i,j}; ysig{i,j} = r * sin(t) + ysig{i,j};
        c1(i,j) = patch(xsig{i,j}, ysig{i,j}, zeros(1,size(xsig{i,j},2)),...
            'linewidth',1.5,'facecolor','none','edgecolor','r');
        xnoise{i,j} = r * cos(t) + xnoise{i,j}; ynoise{i,j} = r * sin(t) + ynoise{i,j};
        c2(i,j) = patch(xnoise{i,j}, ynoise{i,j}, zeros(1,size(xnoise{i,j},2)),...
            'linewidth',1.5,'facecolor','none','edgecolor','g');
    
        tmpsig = unique([round(xsig{i,j})' round(ysig{i,j})'],'rows');
%         tmpsig(:,2) = size(nii{i,j}.img,2) - tmpsig(:,2);
        tmpnoise = unique([round(xnoise{i,j})' round(ynoise{i,j})'],'rows');
%         tmpnoise(:,2) = size(nii{i,j}.img,2) - tmpnoise(:,2);
    
        tmpimgsnr{i,j} = permute(squeeze(nii{i,j}.img(:,:,idx(i,j),1)),[2 1 3]);
        figure(3); t = tmpimgsnr{i,j};
        t(tmpnoise(:,2),tmpnoise(:,1)) = 0; t(tmpsig(:,2),tmpsig(:,1)) = 0;
        imagesc(t)
        
        % snr
        sigInt = tmpimgsnr{i,j}(tmpsig(:,2),tmpsig(:,1),1);
        noiseInt = tmpimgsnr{i,j}(tmpnoise(:,2),tmpnoise(:,1),1);
        snr(i,j) = mean2(sigInt)/std2(noiseInt);
        
        % tsnr
        tmpimgtsnr{i,j} = permute(squeeze(nii{i,j}.img(:,:,idx(i,j),:)),[2 1 3]);
        sigInt = squeeze(sum(sum(tmpimgtsnr{i,j}(tmpsig(:,2),tmpsig(:,1),:),1),2));
        noiseInt = squeeze(sum(sum(tmpimgtsnr{i,j}(tmpnoise(:,2),tmpnoise(:,1),:),1),2));
        tsnr(i,j) = mean(sigInt)/std(sigInt);
    
    end
end
    










%%
addpath(genpath('D:\Preprocessing\NIfTI_20140122\'))

raw_data_fold = {'D:\fMRI\20191011_CW_ofMRI_Parameter_testing\20191011_102736_BEd_ofMRI_20191011_ofMRI_5_1_5\RECON\'};
scannum = {[14:16,18:21,23];34:44};
%{
1 cm:
14:16 250 um thick, 18:21 300 um thick, 23 250 um thick

2 cm:
34:37 250 um thick, 28:44 300 um thick
%}
