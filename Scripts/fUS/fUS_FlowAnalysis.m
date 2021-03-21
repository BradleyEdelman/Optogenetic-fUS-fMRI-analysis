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

MaxI = max(I_neg,I_pos);
I_neg = 10*log10(I_neg./max(max(MaxI)));
I_pos = 10*log10(I_pos./max(max(MaxI)));

I_neg2 = imresize(I_neg,1);
I_pos2 = imresize(I_pos,1);
for  i = 1:size(I_neg2,1)
    for j =1:size(I_neg2,2)
        
        if I_neg2(i,j) > I_pos2(i,j)
            tmp(i,j) = min(I_neg2(:)) - I_neg2(i,j);
        else
            tmp(i,j) = -min(I_pos2(:)) + I_pos2(i,j);
        end
        
%         if tmp(i,j) < 7.5 && tmp(i,j) > -7.5
%             tmp(i,j) = 0;
%         end
        
    end
end
Red = ([linspace(1,0,600)' linspace(0.15,0,600)' linspace(0,0,600)']);
Blue = ([linspace(0,0,600)' linspace(0,.5,600)' linspace(0,1,600)']);
% figure; imagesc(imresize(tmp,1)); caxis([-30 30])
% colormap([flipud(Red);flipud(Blue)])

% I_neg = 10*log10(I_neg./max(max(I_neg)));
% I_pos = 10*log10(I_pos./max(max(I_pos)));

figure;
subplot(2,2,1); imagesc(pdi); title('I total'); set(gca,'yticklabel',[],'xticklabel',[]); colormap(gca,'gray'); caxis([-30 0])
subplot(2,2,2); imagesc(I_pos); title('I Positive'); set(gca,'yticklabel',[],'xticklabel',[]); colormap(gca,'gray'); caxis([-30 0])
subplot(2,2,3); imagesc(I_neg); title('I Negative'); set(gca,'yticklabel',[],'xticklabel',[]); colormap(gca,'gray'); caxis([-30 0])
subplot(2,2,4); imagesc(tmp); title('I Signed'); set(gca,'yticklabel',[],'xticklabel',[]); colormap(gca,[Red;Blue]); caxis([-30 30])

flow.veloc = veloc;
flow.I_neg = I_neg;
flow.I_pos = I_pos;
flow.signed = tmp;


%%
tmp2 = tmp;
tmp2(tmp2 > -4) = 0;
tmp2 = medfilt2(tmp2,[2 2]);
Red = ([linspace(.75,1,600)' linspace(0,1,600)' linspace(0,1,600)']);
figure; imagesc(-tmp2); colormap(flipud(Red));
caxis([0 15])

tmp2 = tmp;
tmp2(tmp2 < 4) = 0;
tmp2 = medfilt2(tmp2,[2 2]);
Blue = ([linspace(1,0,600)' linspace(1,0,600)' linspace(1,.75,600)']);
figure; imagesc(tmp2); colormap(Blue)
caxis([0 15])
        
figure; imagesc(zeros(2,2));
colormap([Red;Blue])
        
        
        