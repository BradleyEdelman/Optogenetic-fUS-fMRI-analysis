function flow = fUS_FlowAnalysis(iqr,pdi)

nfft = 1024;
Fs = 500;
F = linspace(-Fs/2,Fs/2,nfft);

sz = size(iqr);
psd = cell(sz(1:2));
centerF = zeros(sz(1:2));
veloc = zeros(sz(1:2));
I_neg = zeros(sz(1:2));
I_pos = zeros(sz(1:2));


for i = 1:size(iqr,1)
    for j = 1:size(iqr,2)
        
        tmp = squeeze(iqr(i,j,:));
        psd{i,j} = fftshift(abs(fft(tmp,nfft)));
        
        df = Fs/size(psd{i,j},1);
        f = -Fs/2 + df:df:Fs/2;
        
        centerF(i,j) = sum(f'.* psd{i,j}) / sum(psd{i,j});
        veloc(i,j) = centerF(i,j)*1540/(2*15e6)*1000;
        
        I_neg(i,j) = psd{i,j}(1:nfft/2)' * psd{i,j}(1:nfft/2);
        I_pos(i,j) = psd{i,j}(nfft/2 + 1:end)' * psd{i,j}(nfft/2 + 1:end);
        
    end
end

V = f*1540/(2*15e6)*1000;

% global norm makes better visualization 
MaxI = max(I_neg,I_pos);
I_neg_norm = 10*log10(I_neg./max(max(MaxI)));
I_pos_norm = 10*log10(I_pos./max(max(MaxI)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify global normalization same as raw spectrum
from = zeros(size(I_neg)); to = from;
from2 = from; to2 = to;

from(I_neg > I_pos) = 1;
from2(I_neg_norm > I_pos_norm) = 1;
figure(73);
CAX = [1 1 1; 1 0 0 ];
subplot(2,3,1); imagesc(from);
set(gca,'colormap', CAX); caxis([0 1]); title('from trans: psd')
subplot(2,3,2); imagesc(from2); 
set(gca,'colormap', CAX); caxis([0 1]); title('from trans: global norm')

to(I_pos > I_neg) = 1;
to2(I_pos_norm > I_neg_norm) = 1;
CAX = [1 1 1; 0 0 1];
subplot(2,3,4); imagesc(to);
set(gca,'colormap', CAX); caxis([0 1]); title('to trans: psd')
subplot(2,3,5); imagesc(to2);
set(gca,'colormap', CAX); caxis([0 1]); title('to trans: global norm')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_neg2 = imresize(I_neg_norm,1);
I_pos2 = imresize(I_pos_norm,1);
for  i = 1:size(I_neg2,1)
    for j =1:size(I_neg2,2)
        
        % plot dominant flow direction and keep intensity value
        if I_neg2(i,j) > I_pos2(i,j)
            tmp(i,j) = min(I_neg2(:)) - I_neg2(i,j);
        else
            tmp(i,j) = -min(I_pos2(:)) + I_pos2(i,j);
        end
        
    end
end
Red = ([linspace(1,0,600)' linspace(0.15,0,600)' linspace(0,0,600)']);
Blue = ([linspace(0,0,600)' linspace(0,.5,600)' linspace(0,1,600)']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify again signed image contains same "dominant" flow direction
from3 = tmp;
from3(from3 > 0) = 0; from3(from3 < 0) = 1;
figure(73); CAX = [1 1 1; 1 0 0 ];
subplot(2,3,3); imagesc(from3);
set(gca,'colormap', CAX); caxis([0 1]); title('from trans: signed')

to3 = tmp;
to3(to3 < 0) = 0; to3(to3 > 0) = 1;
CAX = [1 1 1; 0 0 1];
subplot(2,3,6); imagesc(to3);
set(gca,'colormap', CAX); caxis([0 1]); title('to trans: signed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2,2,1); imagesc(pdi); title('I total'); set(gca,'yticklabel',[],'xticklabel',[]); colormap(gca,'gray'); caxis([-30 0])
subplot(2,2,2); imagesc(I_pos_norm); title('I Positive'); set(gca,'yticklabel',[],'xticklabel',[]); colormap(gca,'gray'); caxis([-30 0])
subplot(2,2,3); imagesc(I_neg_norm); title('I Negative'); set(gca,'yticklabel',[],'xticklabel',[]); colormap(gca,'gray'); caxis([-30 0])
subplot(2,2,4); imagesc(tmp); title('I Signed'); set(gca,'yticklabel',[],'xticklabel',[]); colormap(gca,[Red;Blue]); caxis([-30 30])

flow.veloc = veloc;
flow.I_neg = I_neg;
flow.I_pos = I_pos;
flow.signed = tmp;
        
        
        