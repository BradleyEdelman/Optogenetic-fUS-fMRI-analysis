
cd('D:\Data\')
addpath(genpath('C:\Users\verasonics\Documents\Vantage-3.4.2-1806191000\Sequences\'))
%%
close all
clear all
folder = 'D:\Data\';

fileA ={'IQdata_11-October-2018_17-32-40.mat'
    };


for i = 1:size(fileA,1)
    fileAtmp = [folder fileA{i}];
    [pdiA{i},timesA{i},blockA{i},psdA{i},velocA{i}] = Clutterfilt(fileAtmp,'filtcutoff',0.3,'velocity',0);
    pdiA{i}=10*log10(pdiA{i}./max(pdiA{i}(:)));
end
figure; imagesc(pdiA{1}); colormap gray; caxis([-35 0])


% L = size(velocA{1}.f,2);
% Fs = 500;
% cutoff = 0.3;
% figure; hold on;
% idx = [10,140];
% x1 = find(velocA{1}.f < -(cutoff*(Fs/2)-.5) & velocA{1}.f > -(cutoff*(Fs/2)+.5)); x1 = velocA{1}.f(x1(end));
% x2 = find(velocA{1}.f > cutoff*(Fs/2)-.5 & velocA{1}.f <cutoff*(Fs/2)+.5); x2 = velocA{1}.f(x2(1));
% maxData = max(velocA{1}.pdipsd{idx(1),idx(2),1});
% patch([x1; x2; x2; x1],[0 0 maxData/2 maxData/2],[.65 .65 .65],'facealpha',1,'linestyle','none')
% area(velocA{1}.f(1:L/2),velocA{1}.pdipsd{idx(1),idx(2),1}(1:L/2),'facecolor','b','edgecolor','b','facealpha',.85)
% area(velocA{1}.f(L/2+1:end),velocA{1}.pdipsd{idx(1),idx(2),1}(L/2+1:end),'facecolor','r','edgecolor','r','facealpha',.85)

fileF = {'IQdata_08-October-2018_16-46-10_Func.mat'
    };
    
for i = 1:size(fileF,1)
    fileFtmp = [folder fileF{i}];
    [pdiF{i},timesF{i},blockF{i},psdF{i},velocF{i}] = Clutterfilt(fileFtmp,'filtcutoff',0.3);
    pdiF{i}(:,:,1:10)=[];

    tmp = findstr('2018',fileF{i});
    time{i} = fileF{i}(tmp+5:tmp+12);
end

%%
remove={[206],[],[],[],[]};
for i = 1:size(pdiF,2)
    blockF{i}(remove{i}) = [];
    pdiF{i}(:,:,remove{i}) = [];
end
%%
for k = 1:size(pdiF,2)
    for j = 1:size(pdiF{k},3)
        pdiF{k}(:,:,j) = imgaussfilt(pdiF{k}(:,:,j),1);
    end
end
%%
% Last index a concatenation of all previous
% sets = size(pdiF,2);
pdiF{sets+1} = cat(3,pdiF{1:sets});
blockF{sets+1} = cat(2,blockF{1:sets});
pdiA{sets+1} = mean(cat(3,pdiA{1:sets}),3);
time{sets+1} = 'Concat';

clear C V
C = cell(1,size(pdiF,2)); V = cell(1,size(pdiF,2));
thresh = cell(1,size(pdiF,2)); cthresh = cell(1,size(pdiF,2));
pdiA2 = cell(1,size(pdiF,2));
for k = 1:size(pdiF,2)
    for i=1:size(pdiF{k},1)
        for j=1:size(pdiF{k},2)

            tmp = squeeze(pdiF{k}(i,j,:));
%             [B,BINT,R,RINT,STATS] = regress(blockF{k}',[1e11*ones(size(pdiF{k},3),1) tmp(:)]);
%             C{k}(i,j) = STATS(1);

            tmp = detrend(squeeze(pdiF{k}(i,j,:)),'linear');
% %             
%             [B,A] = butter(4,0.01,'high');
%             tmp = filtfilt(B,A,tmp);
            
            [r,p] = corrcoef(blockF{k}',tmp(:));
            C{k}(i,j) = r(1,2);
            V{k}(i,j) = var(squeeze(pdiF{k}(i,j,:)));
            
        end
    end
    thresh{k} = std2(C{k})*2.5;
    cthresh{k} = C{k}; 
    cthresh{k}(cthresh{k} < thresh{k} & cthresh{k} > -thresh{k}) = 0;
    % Plot raw and sig (thresholded) correlation map
    figure(100);
    subplot(size(pdiF,2),2,1+2*(k-1)); imagesc(C{k}); title(['Raw Run ' time{k}]); caxis([-1 1])
    subplot(size(pdiF,2),2,2+2*(k-1)); imagesc(cthresh{k}); title(['Thresh Run ' time{k} ': ' num2str(thresh{k})]); caxis([-1 1])
    colormap jet; 
    % Overlay sig correlation map on anatomical
    figure(101);
    pdiA2{k}=pdiA{end}(1:2:end,1:2:end)-1; pdiA2=pdiA;
    pdiA2{k}(cthresh{k} ~= 0) = cthresh{k}(cthresh{k} ~= 0);
    subplot(size(pdiF,2),1,k); imagesc(pdiA2{k})
    title(['R-value: Run ' time{k}])
    caxis([-35 1]);CAX = caxis; Cbar=[gray(abs(CAX(1)-1)*100);jet(abs(CAX(2)+1)*100)]; colormap(Cbar); 

end
figure(100); suptitle('R-values');
% Average of scans
figure(102);
subplot(1,2,1); imagesc(mean(cat(3,C{1:end-1}),3)); caxis([-1 1])
subplot(1,2,2); imagesc(mean(cat(3,cthresh{1:end-1}),3));caxis([-1 1])
colormap jet; suptitle('R-value: Run Average')
figure(103);
C2 = mean(cat(3,C{1:end-1}),3);
thresh2 = std2(C2)*2.5;
cthresh2 = C2; 
cthresh2(cthresh2 < thresh2) = 0;
pdiA22 = pdiA{end}(1:2:end,1:2:end); pdiA22 = pdiA{end};
pdiA22(cthresh2 ~= 0) = cthresh2(cthresh2 ~= 0);
imagesc(pdiA22)
caxis([-35 1]);CAX = caxis; Cbar=[gray(abs(CAX(1)-1)*100);jet(abs(CAX(2)+1)*100)]; colormap(Cbar); 
%%
Y=10:35;
X=[50:60];
base = 1:30;
pditc = cell(1,size(cthresh,2)); start = cell(1,size(cthresh,2));
for k = 1:size(X,1)
    
    figure(104); clf
    for i = 1:size(cthresh,2)
        CORR2tmp = cthresh{i}(Y,X(k,:));
        CORR2tmp = CORR2tmp(:);
        ROI = find(CORR2tmp ~= 0); %ROI = [];
        pditmp = pdiF{i}(Y,X(k,:),:);
        pditmp = reshape(pditmp,size(pditmp,1)*size(pditmp,2),[]);
        
        if isempty(ROI)
            pditc{i} = mean(pditmp(:,:),1);
        else
            pditc{i} = mean(pditmp(ROI,:),1);
        end
        pditc{i} = (pditc{i}-mean(pditc{i}(base)))/(mean(pditc{i}(base)))*100;
        pditc{i} = detrend(pditc{i},'linear');
        subplot(size(cthresh,2),1,i); plot(pditc{i},'r','linewidth',2); hold on; plot(25*blockF{i},'k','linewidth',1)
        set(gca,'ylim',[-25 65]); grid minor; title(['Timecourse: Run ' time{i}]) 

        start{i} = find(diff(blockF{i}) == 1)+1;
        
    end
%     kk = 1; file = [folder 'Timecourse_X_' num2str(X(k,1)) '-' num2str(X(k,end)) '_Y_' num2str(Y(k,1)) '-' num2str(Y(k,end)) '_' num2str(kk) '.fig'];
%     while isequal(exist(file,'file'),2); kk = kk +1; file = [folder 'Timecourse_X_' num2str(X(k,1)) '-' num2str(X(k,end)) '_Y_' num2str(Y(k,1)) '-' num2str(Y(k,end)) '_' num2str(kk) '.fig']; end
%     savefig(figure(104),file)
    
    figure(105); clf
    trialstd = cell(1,size(start,2)); trial2 = cell(1,size(start,2));
    for i = 1:size(start,2)
        
        trial = cell(1,size(start{i},2));
        for j = 1:size(start{i},2)
            trial{i}(1,:,j) = pditc{i}(start{i}(j)-11:start{i}(j)+19);
        end
        trialstd{i} = std(trial{i},0,3);
        trial2{i} = mean(trial{i},3);

        subplot(size(start,2),1,i); hold on; grid on; title(['Timecourse: Run ' time{i}]) 
        errorbar(-11:19,trial2{i},trialstd{i},'k','linewidth',1.5)
        minData = min(trial2{i});
        maxData = max(trial2{i});
        patch([0; 10; 10; 0],[minData-10 minData-10 minData-7.5 minData-7.5],...
            'b','facealpha',0.75,'linestyle','none')
        axis([-10 20 minData-10 maxData+20])
    end
%     kk = 1; file = [folder 'Evoked_X_' num2str(X(k,1)) '-' num2str(X(k,end)) '_Y_' num2str(Y(k,1)) '-' num2str(Y(k,end)) '_' num2str(kk) '.fig'];
%     while isequal(exist(file,'file'),2); kk = kk +1; file = [folder 'Evoked_X_' num2str(X(k,1)) '-' num2str(X(k,end)) '_Y_' num2str(Y(k,1)) '-' num2str(Y(k,end)) '_' num2str(kk) '.fig']; end
%     savefig(figure(105),file)
    
end
%%
% Overlay max and average response during block stimulation
base = 1:30;
for k = 1:size(pdiF,2)
    
    for i = 1:size(pdiF{k},1)
        for j = 1:size(pdiF{k},2)
    
            tmp = squeeze(pdiF{k}(i,j,:));
            tmp = (tmp-mean(tmp(base)))/(mean(tmp(base)))*100;
            tmp = detrend(tmp,'linear');
            start = find(diff(blockF{k}) == 1)+1;
            for l = 1:size(start,2)
                trial = tmp(start(l)-11:start(l)+19);
                maxtrial(l) = max(trial(11:20));
                meantrial(l) = mean(trial(11:20));
            end
            meanmax{k}(i,j) = mean(maxtrial);
            meanmean{k}(i,j) = mean(meantrial);
        end
    end
    pdiA3{k}=pdiA{end}(1:2:end,1:2:end)-1; pdiA3{k} = pdiA{end}
    pdiA3{k}(cthresh{k} ~= 0) = meanmax{k}(cthresh{k} ~= 0);
    figure(106)
    subplot(size(pdiF,2),1,k); imagesc(pdiA3{k})
    title(['Max Response: Run ' time{k}])
    caxis([-26 30]);CAX = caxis; Cbar=[gray(abs(CAX(1)-1)*100);jet(abs(CAX(2)+1)*100)]; colormap(Cbar);
    pdiA3{k}(cthresh{k} ~= 0) = meanmean{k}(cthresh{k} ~= 0);
    figure(107)
    subplot(size(pdiF,2),1,k); imagesc(pdiA3{k})
    title(['Mean Response: Run ' time{k}])
    caxis([-26 10]);CAX = caxis; Cbar=[gray(abs(CAX(1)-1)*100);jet(abs(CAX(2)+1)*100)]; colormap(Cbar);
end

            
            
            

%%
k = 1; file = [folder 'R-value_' num2str(k) '.fig'];
while isequal(exist(file,'file'),2); k = k +1; file = ['R-values ' num2str(k) '.fig']; end
savefig(figure(100),file)

k = 1; file = [folder 'R-value_ANAT_' num2str(k) '.fig'];
while isequal(exist(file,'file'),2); k = k +1; file = [folder 'R-value_ANAT_' num2str(k) '.fig']; end
savefig(figure(101),file)

k = 1; file = [folder 'R-value_AVE_' num2str(k) '.fig'];
while isequal(exist(file,'file'),2); k = k +1; file = [folder 'R-value_AVE_' num2str(k) '.fig']; end
savefig(figure(102),file)

k = 1; file = [folder 'R-values_ANAT_AVE_' num2str(k) '.fig'];
while isequal(exist(file,'file'),2); k = k +1; file = [folder 'R-values_ANAT_AVE_' num2str(k) '.fig']; end
savefig(figure(103),file)





%%
VOX =10*[.1 .1 .1];
PDI = permute(PDI, [1 2 4 3]);
nii = make_nii(PDI,VOX,[0 0 0]);
        save_nii(nii, 'PDIFunc2.nii');
% %         spm_hwrite('PDI3',[118 130 1 180],VOX,1,4,0,[0 0 0],'spm compatible');






