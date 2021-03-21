D.IQR=IQSAVE;
D.times=times;

D.block=[zeros(1,30),ones(1,10),zeros(1,20),ones(1,10),zeros(1,20),ones(1,10),zeros(1,20),ones(1,10),zeros(1,20)];
D.notes = '';


%%
clear PDI
for k=1:size(D.IQR,4)
    k
    IQR=D.IQR(:,:,:,k);
    
%     if k<10
%         eval(['IQR = frame_00' num2str(k) ';']);
%     else
%         eval(['IQR = frame_0' num2str(k) ';']);
%     end
    
    for i=1:size(IQR,1)
        for j=1:size(IQR,2)
            IQR(i,j,:)=IQR(i,j,:)-IQR(i,j,1); 
        end
    end
    IQR1=IQR;
%     [nz,nx,nt]=size(IQR);
% 
%     m=zeros(nx*nz,nt);
%     i=0;
%     for ix=1:nx
%         for iz=1:nz
%             i=i+1;
%             m(i,:)=squeeze(IQR(iz,ix,:)-mean(IQR(iz,ix,:)));
%         end
%     end
%     clear IQR
% 
% 
%     [U,S,V]=svd(m,'econ');
%     clear m
%     
%     for l=1:size(V,2)
%         [pxx,f]=periodogram(V(:,l),tukeywin(200,.2),'centered');
%         PSD(:,l,k)=10*log10(pxx);
%     end
%     
%     nfilt=25;
%     Sf=S;
%     Sf(1:nfilt,1:nfilt)=0;
%     Sf(end)=0;
%     mf=U*Sf*V';
%     clear U S Sf V
% 
%     i=0;
%     IQR1=zeros(nz,nx,nt);
%     for ix=1:nx
%         for iz=1:nz
%             i=i+1;
%             IQR1(iz,ix,:)=squeeze(mf(i,:));
%         end
%     end
%     clear mf

    [B,A]=butter(4,0.3,'high');
    IQR2=filter(B,A,IQR1,[],3);
    IQR2=IQR2(:,:,4:end);           % the first 4 temporal samples are eliminates (filter oscilations)
    PDI(:,:,k)=mean(abs(IQR2).^2,3);     % computing the intensity of the blood signal the 
    
end
PDI(:,:,1:10)=[];
% 
% PDIdisplay=zeros(size(PDI,1),size(PDI,2),1,size(PDI,3));
% PDInii=zeros(size(PDI,1),size(PDI,2),1,size(PDI,3));
clear PDIdisplay PDInii
for i=1:size(PDI,3)
    PDItmp=PDI(:,:,i);
    
    PDIdisplay(:,:,1,i)=10*log10(PDItmp./max(PDItmp(:)));
    PDInii(:,:,1,i)=PDItmp;
end

figure; imagesc(PDIdisplay(:,:,1,15));colormap gray; colorbar

block=[zeros(1,30),ones(1,10),zeros(1,20),ones(1,10),zeros(1,20),ones(1,10),zeros(1,20),ones(1,10),zeros(1,20)];
% block = D.block;

remove=[];
block(remove)=[];
PDI(:,:,remove)=[];
% PDI=PDI-repmat(mean(PDI,3),1,1,size(PDI,3));

clear CORR VAR
for i=1:size(PDIdisplay,1)
    for j=1:size(PDI,2)
        
        [B,BINT,R,RINT,STATS] = regress(block',[ones(size(PDI,3),1) squeeze(PDI(i,j,:,:))]);
        CORR(i,j)=STATS(1);
        VAR(i,j)=var(squeeze(PDI(i,j,:)));
    end
end
CORR2=CORR; CORR2(CORR2<0.2)=0;
figure; subplot(1,2,1); imagesc(CORR); subplot(1,2,2); imagesc(CORR2)
colormap jet
%%
X=45;
Y=15;
figure; subplot(1,2,1); plot(block*1e10,'k'); hold on; plot(squeeze(PDI(Y,X,:)),'r'); grid minor
title(['X: ' num2str(X) ' Y: ' num2str(Y)]);
subplot(1,2,2); plot(block*1e10,'k'); hold on; plot(squeeze(mean(mean(PDI,1),2)),'r'); grid minor
title('Global Mean')
%%
VOX=[1 1 1];
nii = make_nii(PDI,VOX,[0 0 0], 4);
        save_nii(nii, 'PDI');
% %         spm_hwrite('PDI3',[118 130 1 180],VOX,1,4,0,[0 0 0],'spm compatible');






