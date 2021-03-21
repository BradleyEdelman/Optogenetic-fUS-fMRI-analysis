nTRs=360;
t=0:0.1:540;
box=[zeros(1,900) repmat([ones(1,120) zeros(1,780)],1,5) 0];
box=[zeros(1,900) repmat([ones(1,120) zeros(1,795)],1,5) 0]; box = box(1:5401);

nTRs = 240;
t = 0:0.1:360;
box = [zeros(1,600) repmat([ones(1,120) zeros(1,465)],1,5) zeros(1,250)]; box = box(1:3601);


figure(250); subplot(1,3,1); plot(t,box);axis([0 540 0 1.5]); xlabel('sec'); title('Stimulus Block Design')

T0=0; n=4; lamda=2;
hrf=((t-T0).^(n-1)).*exp(-(t-T0)/lamda)/((lamda^n)*factorial(n-1));
figure(250); subplot(1,3,2); plot(hrf(1:240));
xlabel('1/10th seconds'); title('HRF')

B=conv(hrf,box)/10;
tp=0:.1:2400;
for i=1:nTRs
    N(i)=B(i*15);
end
figure(250); subplot(1,3,3); plot(N); axis([0 nTRs 0 1.5]); xlabel('Acq Number'); title('Stimulus Block * HRF')
% N = circshift(N,9,2)

clear x
x(:,1) = N';
x(:,2) = ones(nTRs,1);
x(:,3) = linspace(1,nTRs,nTRs)';
x(:,4) = [linspace(1,nTRs,nTRs/2) linspace(nTRs,1,nTRs/2)];
x(:,5) = [linspace(1,nTRs,nTRs/3) linspace(nTRs,1,nTRs/3) linspace(1,nTRs,nTRs/3)];
x(:,6) = repmat([linspace(1,nTRs,nTRs/4) linspace(nTRs,1,nTRs/4)],1,2);

x(:,3:end)=x(:,3:end)/nTRs;
X = [];
for i = 1:size(pdiF,2)
    X = blkdiag(X,x);
end
figure(251);
c=repmat([1;zeros(size(x,2)-1,1)],size(pdiF,2),1);
subplot(2,1,1); 
b=bar(.5:1:size(c,1),c,'k','barwidth',.5); grid minor; title('Design Matrix')
set(gca,'position',[.13 .80 .775 .15],'ylim',[-2 2],'xlim',[0 size(c,1)],...
    'xtick',0:1:12,'xticklabel',[],'yticklabel',[])
subplot(2,1,2);
image(X*64);
set(gca,'position',[.13 .1 .775 .65],'xticklabel',[],'yticklabel',[1 360],'ytick',[1 size(N,2)])
ylabel('Acq Number');
colormap(flipud(gray));
nTR = nTRs * size(pdiF,2);





%%%
tdist2T = @(t,DF) (1-betainc(DF/(DF+t^2),DF/2,0.5)); % 2-tailed t-distribution
tdist1T = @(t,DF) 1-(1-tdist2T(t,DF))/2; % 1-tailed t-distribution

PDI = cat(3,pdiF{:});
clear Beta
for i=1:size(PDI,1)
    for j=1:size(PDI,2)

        Y = squeeze(PDI(i,j,:))';
        
        DF = nTR - 2;
        
        beta_hat=inv(X'*X)*X'*Y';
        Beta(i,j,:) = reshape(beta_hat,1,1,size(c,1));
%         Beta(i,j,2) = beta_hat(2);
%         Beta(i,j,3) = beta_hat(3);
        Var_e=(Y'-X*beta_hat)'*(Y'-X*beta_hat)/DF;
        
        %Hypothesis testing; Compute the t statistic
        t_stat(i,j)=c'*beta_hat/sqrt(Var_e*c'*inv(X'*X)*c);
        
        p2tail(i,j) = 1- tdist2T(t_stat(i,j),DF);
        p1tail(i,j) = 1- tdist1T(t_stat(i,j),DF);

    end
end
figure(252);
for i = 1:size(c,1)
    d1 = round(size(c,1)/size(pdiF,2));
    d2 = round(size(c,1)/d1);
    subplot(d2,d1,i); imagesc(Beta(:,:,i));
    set(gca,'yticklabel',[],'xticklabel',[])
end
colormap jet

Pthresh = 0.005; Maxt = ceil( max(t_stat(:)));
CAX = [-35-Maxt Maxt]; Cbar=[gray(36*50);winter(abs(CAX(2)+1)*50);autumn(abs(CAX(2)+1)*50)];
figure(253);
subplot(3,2,1); imagesc(p2tail); caxis([0 1]); title('2 tail - P values'); colormap(gca,'jet'); caxis([0 1])
subplot(3,2,2); imagesc(p1tail); caxis([0 1]); title('1 tail - P values'); colormap(gca,'jet'); caxis([0 1])

pdiA2 = pdiA{end} - Maxt; pdiA2(p2tail < Pthresh) = t_stat(p2tail < Pthresh); 
subplot(3,2,3); imagesc(pdiA2); colormap(gca,Cbar); caxis(CAX); title(['2 tail - P < ' num2str(Pthresh)])
pdiA2 = pdiA{end} - Maxt; pdiA2(p1tail < Pthresh) = t_stat(p1tail < Pthresh);
subplot(3,2,4); imagesc(pdiA2); colormap(gca,Cbar); caxis(CAX); title(['1 tail - P < ' num2str(Pthresh)])

[c_pvalues, c_alpha, h] = fwer_sidak(reshape(p2tail,size(p2tail,1)*size(p2tail,2), []), 0.05);
c_p2values = reshape(c_pvalues,size(p2tail,1),size(p2tail,2));
pdiA2 = pdiA{end} - Maxt; pdiA2(c_p2values < Pthresh) = t_stat(c_p2values < Pthresh); 
subplot(3,2,5); imagesc(pdiA2); colormap(gca,Cbar); caxis(CAX); title(['2 tail - P < ' num2str(Pthresh) ',FWE Corr'])

[c_pvalues, c_alpha, h] = fwer_sidak(reshape(p1tail,size(p1tail,1)*size(p1tail,2), []), 0.05);
c_p1values = reshape(c_pvalues,size(p1tail,1),size(p1tail,2));
pdiA2 = pdiA{end} - Maxt; pdiA2(c_p1values < Pthresh) = t_stat(c_p1values < Pthresh); 
subplot(3,2,6); imagesc(pdiA2); colormap(gca,Cbar); caxis(CAX); title(['1 tail - P < ' num2str(Pthresh) ',FWE Corr'])



sets = 1;
blockF{1} = N;




