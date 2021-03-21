function [pdi,times,block,PSD,velocity] = Clutterfilt(file,varargin)

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
    iqr = Data.(fieldNames{end});
    pdi = zeros(size(iqr,1),size(iqr,2),1);
    times = [];
    block = [];
   
elseif isfield(Data,'D') % low resolution functional
    
    iqr = Data.D.IQR;
    pdi = zeros(size(iqr,1),size(iqr,2),size(iqr,4));
    times = Data.D.times;
    block = Data.D.block;
    
end


PSD = [];
iqrv = zeros(size(iqr,1),size(iqr,2),size(iqr,3),size(iqr,4));
for k=1:size(iqr,4)
k
    iqrtmp=iqr(:,:,:,k);
    
    for i=1:size(iqrtmp,1)
        for j=1:size(iqrtmp,2)
            iqrtmp(i,j,:)=iqrtmp(i,j,:)-iqrtmp(i,j,1); 
        end
    end
    iqr1=iqrtmp;
    
    if strcmp(p.method,'svd')
    
        [nz,nx,nt]=size(iqr);

        m=zeros(nx*nz,nt);
        i=0;
        for ix=1:nx
            for iz=1:nz
                i=i+1;
                m(i,:)=squeeze(iqr(iz,ix,:)-mean(iqr(iz,ix,:)));
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
        iqr1=zeros(nz,nx,nt);
        for ix=1:nx
            for iz=1:nz
                i=i+1;
                iqr1(iz,ix,:)=squeeze(mf(i,:));
            end
        end
        clear mf
        
    end
    
    [B,A] = butter(4,0.3,'high');
    iqr2 = filter(B,A,iqr1,[],3);
    iqrv(:,:,:,k) = iqr2;
    iqr2 = iqr2(:,:,4:end);                % the first 4 temporal samples are eliminates (filter oscilations)
    pdi(:,:,k) = mean(abs(iqr2).^2,3);     % computing the intensity of the blood signal
    
end


if isfield(p,'velocity') && isequal(p.velocity,1)
    
    nfft = 1024;
    
    pdiv = zeros(size(pdi,1),size(pdi,2),size(pdi,3));
    centerF = zeros(size(pdi,1),size(pdi,2),size(pdi,3));
    veloc = zeros(size(pdi,1),size(pdi,2),size(pdi,3));
    pdipsd = cell(size(pdi,1),size(pdi,2),size(pdi,3));
    pdiplot = zeros(size(pdi,1),size(pdi,2),size(pdi,3));
    
    for k = 1:size(pdi,3)
        
        pdiv(:,:,k) = pdi(:,:,k);
        pdiv(:,:,k) = 10*log10(pdiv(:,:,k)./max(max(pdiv(:,:,k))));
    
        for i = 1:size(pdiv,1)
            for j = 1:size(pdiv,2)
                
                if pdiv(i,j,k) > -25 % Threshold PDI power at -25 dB
                    
                    tmp = squeeze(iqrv(i,j,:));
                    pdipsd{i,j,k} = fftshift(abs(fft(tmp,nfft)));
                    
                    df = 500/size(pdipsd{i,j,k},1);
                    f = -500/2+df:df:500/2;
                    
                    centerF(i,j,k) = sum(f'.* pdipsd{i,j,k}) / sum(pdipsd{i,j,k});
                    veloc(i,j,k) = centerF(i,j,k)*1540/(2*15e6)*1000;
                    
                else
                    
                    pdipsd{i,j,k} = zeros(nfft,1);
                    centerF(i,j,k) = 0;
                    veloc(i,j,k) = 0;
                    
                end
                
            end
        end
     
        figure(10); clf
        imagesc(veloc(:,:,k));
        cmap=jet(256); cmap(90:158,:)=repmat([0 0 0],69,1); colormap(cmap)
    
        minData = min(min(veloc));
        pditmp = pdiv(:,:,k) + minData;
        pditmp(veloc < -1.5 | veloc > 1.5) = veloc(veloc < -1.5 | veloc > 1.5);
        pdiplot(:,:,k) = pditmp;
        figure(11); clf; imagesc(pdiplot)
        caxis([-25+minData max(max(veloc))]); CAX = caxis;
        Cbar=[gray(round(abs(CAX(1)-minData)*100));jet(round(abs(CAX(2)-minData)*100))];
        colormap(Cbar);
        
        [filepath,name,ext] = fileparts(file);
        kk = 1; file = [filepath '/Veloc_' num2str(kk) '.fig'];
        while isequal(exist(file,'file'),2); kk = kk +1; file = [filepath '/Veloc_' num2str(kk) '.fig']; end
%         savefig(figure(11),file)
        
    end
    
    velocity.pdipsd = pdipsd;
    velocity.centerF = centerF;
    velocity.veloc = veloc;
    velocity.pdiv = pdiv;
    velocity.pdiplot = pdiplot;
    velocity.f = f;
    
else
    
    velocity = [];
    
end








