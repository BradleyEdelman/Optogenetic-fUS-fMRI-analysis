function [SNR,tSNR,fusroi] = fusi_compute_snr(Data)

% Data = squeeze(Data(:,:,10,:));
Anat = 10*log10(Data(:,:,1)./(max(max(Data(:,:,1)))));
f=figure; imagesc(abs(Anat)); axis tight; colormap gray

[xnoise, ynoise] = getpts(f);
t = linspace(0, 2*pi); r = 2;
xnoise = r * cos(t) + xnoise; ynoise = r * sin(t) + ynoise;
c2 = patch(xnoise, ynoise, zeros(1,size(xnoise,2)),...
    'linewidth',1.5,'facecolor','none','edgecolor','g');
tmpnoise = unique([round(xnoise)' round(ynoise)'],'rows');
tmpnoise((sum(tmpnoise<=0,2)>0),:) = []; 


% Load segmentation
roifolder = 'D:\Data_Processed\segment_ofUS\fus_20191205\';
roi = dir(roifolder);
roi(contains({roi.name},'.fig')) = [];
roi(~contains({roi.name},'R')) = [];
test = zeros(size(Data(:,:,1)));
for i_roi = 1:size(roi,1)
    roifile = [roifolder roi(i_roi).name];
    load(roifile)
    cdata = imresize(ROI.cdata,.5);
    pts{i_roi} = find(cdata);
    test(pts{i_roi}) = i_roi;
end
figure(100); imagesc(test)

% [xsig, ysig] = getpts(f);
% 
% xsig = r * cos(t) + xsig; ysig = r * sin(t) + ysig;
% c1 = patch(xsig, ysig, zeros(1,size(xsig,2)),...
%     'linewidth',1.5,'facecolor','none','edgecolor','r');
% tmpsig = unique([round(xsig)' round(ysig)'],'rows');
% tmpsig((sum(tmpsig<=0,2)>0),:) = [];

[b,a] = butter(5,0.01/((1/1.5)/2),'high');


Data2 = reshape(Data,size(Data,1)*size(Data,2),1,size(Data,3));
for i_roi = 1:size(roi,1)

    roi_data = Data2(pts{i_roi},:,:);
    noise_data = Data(tmpnoise(:,2),tmpnoise(:,1),:);
%     figure(50); subplot(5,4,i_roi); plot(detrend(squeeze(mean(roi_data,1)),'linear'));
%     figure(50); subplot(5,4,i_roi); plot(squeeze(mean(roi_data,1)));

    roi_data2 = filtfilt(b,a,squeeze(mean(roi_data,1)));
    noise_data2 = filtfilt(b,a,squeeze(mean(mean(noise_data,1),2)));

    figure(50); subplot(5,4,i_roi); cla
    plot(roi_data2,'b'); hold on
    plot(noise_data2,'r');
    drawnow
    % snr (first image)
    S = mean(roi_data(:,:,1));
    N = std2(noise_data(:,:,1));
    SNR(i_roi) = S/N;
    
    % tsnr (all images)
    S = mean(squeeze(mean(roi_data(:,:,1:30),1)));
    N = std(squeeze(mean(roi_data(:,:,1:30),1)));
%     N = std(squeeze((mean(mean(noise_data,1),2))));
    tSNR(i_roi) = S/N;
end

fusroi = {'LCPu' 'LS1Dz' 'LS1FL' 'LS1HL' 'LM1' 'LS2' 'LM2' 'LCg1' 'LCg2'...
    'RCg1' 'RCg2' 'RM2' 'LS1BF' 'RM1' 'RCPu' 'RS1HL' 'RS1FL' 'RS1Dz' 'RS1BF' 'RS2'};
