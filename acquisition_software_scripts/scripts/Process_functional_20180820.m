
cd('/oak/stanford/groups/ljinhy/bjedelma')
addpath(genpath('/oak/stanford/groups/ljinhy/bjedelma/US_IMAGING/scripts'))
%%
% fileAnat = '/home/bradley/oak/bjedelma/US_IMAGING/IQdata_14-August-2018_16-20-09.mat';
fileAnat = '/oak/stanford/groups/ljinhy/bjedelma/US_IMAGING/IQdata_13-August-2018_11-11-58.mat';
[pdiAnat,timesAnat,blockAnat,psdAnat] = Clutterfilt(fileAnat,'filtcutoff',0.3);
pdiAnat=10*log10(pdiAnat./max(pdiAnat(:)));

fileF = '/oak/stanford/groups/ljinhy/bjedelma/US_IMAGING/IQdata_13-August-2018_11-11-58_RHL_STIM_25uA_GOOD_THAL.mat';
[pdiF,timesF,blockF,psdF] = Clutterfilt(fileF,'filtcutoff',0.3);
pdiF(:,:,1:10)=[];

remove=[41,62,88];
blockF(remove)=[];
pdiF(:,:,remove)=[];

clear C V
for i=1:size(pdiF,1)
    for j=1:size(pdiF,2)
        [B,BINT,R,RINT,STATS] = regress(blockF',[ones(size(pdiF,3),1) squeeze(pdiF(i,j,:,:))]);
        C(i,j)=STATS(1);
        V(i,j)=var(squeeze(pdiF(i,j,:)));
    end
end
thresh=std2(C)*2;
cthresh=C; 
cthresh(cthresh<thresh)=0;
figure; subplot(1,2,1); imagesc(C); subplot(1,2,2); imagesc(cthresh)
colormap jet

% Superimpose on anat
pdiAnat2=pdiAnat(1:2:end,1:2:end);
pdiAnat2(cthresh ~= 0) = cthresh(cthresh ~= 0);
figure; imagesc(pdiAnat2)
caxis([-25 .5]);CAX = caxis; Cbar=[gray(abs(CAX(1))*100);jet(CAX(2)*100)]; colormap(Cbar); 

CORR2tmp = cthresh(1:20,40:50);
CORR2tmp = CORR2tmp(:);
ROI = find(CORR2tmp ~= 0);
pditmp = pdiF(1:20,40:50,:);
pditmp = reshape(pditmp,size(pditmp,1)*size(pditmp,2),[]);

pditc = mean(pditmp(ROI,:),1);
pditc = (pditc-mean(pditc(1:20)))/(mean(pditc(1:20)))*100;
pditc = detrend(pditc,'linear');
figure; plot(pditc,'r','linewidth',2); hold on; plot(25*blockF,'k','linewidth',1)
set(gca,'ylim',[-50 100])

%%
X=70;
Y=50;
figure; subplot(1,2,1); plot(block*1e10,'k'); hold on; plot(squeeze(PDI(Y,X,:)),'r'); grid minor
title(['X: ' num2str(X) ' Y: ' num2str(Y)]);
subplot(1,2,2); plot(block*1e10,'k'); hold on; plot(squeeze(mean(mean(PDI,1),2)),'r'); grid minor
title('Global Mean')
%%
VOX =10*[.1 .1 .1];
pdiF = permute(pdiF, [1 2 4 3]);
nii = make_nii(pdiF,VOX,[0 0 0]);
        save_nii(nii, 'PDIFunc.nii');
% %         spm_hwrite('PDI3',[118 130 1 180],VOX,1,4,0,[0 0 0],'spm compatible');






