function [Data,PDI,PSD,IQR2] = Clutterfilt_preproc(file,varargin)

if ~rem(size(varargin,2),2) == 0
    error('Odd number of inputs/n');
else
    for i = 1:2:size(varargin,2)
        p.(varargin{i}) = varargin{i+1};
    end
end

if ~isfield(p,'filtcutoff')
    error('need to define filter cutoff/n');
elseif p.filtcutoff > 1
    p.filtcutoff = p.filtcutoff/(500/2);
end

if ~isfield(p,'method')
    p.method = 'nosvd';
else
    if strcmp(p.method,'svd')
        if ~isfield(p,'svdcutoff')
            error('Need to specify cutoff value for SVD method/n');
        end
    end
end

Data = load(file);
if isfield(Data,'frame_001') % High resolution anatomical
    
    fieldNames = fieldnames(Data);
    IQR = Data.(fieldNames{end});
    PDI = zeros(size(IQR,1),size(IQR,2),1);
    
%     for  i = 1:size(fieldNames,1)
%         IQR(:,:,:,i) = Data.(fieldNames{i});
%     end
%     clear Data
%     p.method = 'nosvd';
%    
elseif isfield(Data,'Data') % low resolution functional
    
    Data = Data.Data;
    IQR = Data.D.IQR;
    PDI = zeros(size(IQR,1),size(IQR,2),size(IQR,4));
    
end


PSD = [];
for k=1:size(IQR,4)
    k
    IQRtmp=IQR(:,:,:,k);
    
    for i=1:size(IQRtmp,1)
        for j=1:size(IQRtmp,2)
            IQRtmp(i,j,:)=IQRtmp(i,j,:)-IQRtmp(i,j,1); 
        end
    end
    IQR1=IQRtmp;
    
    if strcmp(p.method,'svd')
    
        [nz,nx,nt]=size(IQR);

        m=zeros(nx*nz,nt);
        i=0;
        for ix=1:nx
            for iz=1:nz
                i=i+1;
                m(i,:)=squeeze(IQR(iz,ix,:)-mean(IQR(iz,ix,:)));
            end
        end
        clear IQR

        [U,S,V]=svd(m,'econ');
        clear m

        for l=1:size(V,2)
            [pxx,f]=periodogram(V(:,l),tukeywin(200,.2),'centered');
            PSD(:,l,k)=10*log10(pxx);
        end

        nfilt=p.svdcutoff;
        Sf=S;
        Sf(1:nfilt,1:nfilt)=0;
        Sf(end)=0;
        mf=U*Sf*V';
        clear U S Sf V

        i=0;
        IQR1=zeros(nz,nx,nt);
        for ix=1:nx
            for iz=1:nz
                i=i+1;
                IQR1(iz,ix,:)=squeeze(mf(i,:));
            end
        end
        clear mf
        
    end
    
    [B,A]=butter(4,0.3,'high');
    IQR2=filter(B,A,IQR1,[],3);
%     IQR2=IQR2(:,:,4:end);           % the first 4 temporal samples are eliminates (filter oscilations)
    PDI(:,:,k)=mean(abs(IQR2(:,:,4:end)).^2,3);     % computing the intensity of the blood signal the 
    
end

% 
% figure;
% for  i = 1:size(PDI,3)
%     pdi2(:,:,i) = 10*log10(PDI(:,:,i)./max(max(PDI(:,:,i))));
%     subplot(ceil(size(PDI,3)/5),5,i);
%     imagesc(pdi2(:,:,i))
%     caxis([-35 0])
%     colormap gray
%     drawnow
% end
