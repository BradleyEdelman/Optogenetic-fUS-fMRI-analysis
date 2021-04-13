function GLM = GLM_fMRI(func_fold,anat_fold,srfile,param)


if isunix; slash = '/'; elseif ispc; slash = '\'; end

dummy = param.Dummy;
order = param.order;
template = param.template;

nii = load_untouch_nii(srfile);

if ~isempty(anat_fold)
    
    
    ave_anat = [anat_fold 'ave_anatomy.nii'];
    if exist(ave_anat,'file') && isequal(template,1)
        
        anii = load_untouch_nii(ave_anat);
        
        ave_anat_mask = [anat_fold 'ave_anatomy_mask.nii'];
        if exist(ave_anat_mask,'file')
            mnii = load_untouch_nii(ave_anat_mask);
            mask = mnii.img;
        else
            mask = ones(size(anii.img));
        end
        
    else
    
        anat = [anat_fold anat_fold(end-3:end-2) '_anatomy.nii'];
        anii = load_untouch_nii(anat);
        anii.img = flipud(permute(anii.img,[2 1 3]));
        
        anat_mask = [anat_fold anat_fold(end-3:end-2) '_anatomy_mask.nii'];
        if exist(anat_mask,'file')
            mnii = load_untouch_nii(anat_mask);
            mnii.img = flipud(permute(mnii.img,[2 1 3]));
            mask = mnii.img;
        else
            mask = ones(size(anii.img));
        end
        
    end
    
else
    
    anii = nii;
    anii.img = anii.img(:,:,:,1);
    mask = ones(size(anii.img,1),size(anii.img,2),size(anii.img,3));
    
end

Dir = pwd;
if isunix
    load([Dir(1:41) 'Preprocessing/4th_order_gam_Design_fMRI.mat'])
elseif ispc
    load([Dir(1:3) 'Preprocessing\4th_order_gam_Design_fMRI.mat'])
end
switch order
    case 1
        X(2:4,:) = [];
    case 2
        X(3:4,:) = [];
    case 3
        X(4,:) = [];
    case 4
end

X = X(dummy+1:end,:);
nTRs = size(nii.img,4);

% Define t-ditributions
tdist2T = @(t,DF) (1-betainc(DF/(DF+t^2),DF/2,0.5)); % 2-tailed
tdist1T = @(t,DF) 1-(1-tdist2T(t,DF))/2; % 1-tailed

c = [0;1;zeros(size(X,2)-2,1)]; % utilize 2nd order gamma
beta = zeros(size(nii.img,1),size(nii.img,2),size(nii.img,3),size(c,1));
t_stat = zeros(size(nii.img,1),size(nii.img,2),size(nii.img,3));
p2tail = zeros(size(nii.img,1),size(nii.img,2),size(nii.img,3));
p1tail = zeros(size(nii.img,1),size(nii.img,2),size(nii.img,3));
for i = 1:size(nii.img,1)
    for j = 1:size(nii.img,2)
        for k = 1:size(nii.img,3)

            Y = squeeze(nii.img(i,j,k,:))';

            DF = nTRs - 2;

            beta_hat=inv(X'*X)*X'*Y';
            beta(i,j,k,:) = reshape(beta_hat,1,1,size(c,1));
            Var_e=(Y'-X*beta_hat)'*(Y'-X*beta_hat)/DF;

            %Hypothesis testing; Compute the t statistic
            t_stat(i,j,k)=c'*beta_hat/sqrt(Var_e*c'*inv(X'*X)*c);
            
            if isnan(t_stat(i,j,k))
                p2tail(i,j,k) = 1;
                p1tail(i,j,k) = 1;
            else
                p2tail(i,j,k) = 1 - tdist2T(t_stat(i,j,k),DF);
                p1tail(i,j,k) = 1 - tdist1T(t_stat(i,j,k),DF);
            end
            
        end
    end
end

GLM.X = X;
GLM.beta = beta;
GLM.t_stat = t_stat;
GLM.p1tail = p1tail;
GLM.p2tail = p2tail;

h1 = figure(10); subplot(2,1,1); 
b=bar(.5:1:size(c,1),c,'k','barwidth',.5); grid minor; title('Design Matrix')
set(gca,'position',[.13 .80 .775 .15],'ylim',[-2 2],'xlim',[0 size(c,1)],...
    'xtick',0:1:12,'xticklabel',[],'yticklabel',[])
subplot(2,1,2); image(X*64);
set(gca,'position',[.13 .1 .775 .65],'xticklabel',[],'yticklabel',[1 size(X,1)],'ytick',[1 size(X,2)])
ylabel('Acq Number');
colormap gray

Pthresh = 0.005; maxt = 20;

plotmin = -100; plotmax = 0;
Cbar = [gray(abs(plotmax - plotmin)*100); fireice(2*maxt*100)]; CAX = [plotmin - maxt maxt];
Cbarp = [jet(round(1.1*100));[0 0 0]];
h2 = figure(11); h3 = figure(12); h4 = figure(13); h5 = figure(14);

% Resize time-series param to anat sampling
p2tmp = imresize3(p2tail,size(anii.img));
p2tmp(mask == 0) = 1.1;
p1tmp = imresize3(p1tail,size(anii.img));
p1tmp(mask == 0) = 1.1;
ttmp = imresize3(t_stat,size(anii.img));
ttmp(mask == 0) = 0;

% rescale dynamic range of anat for vis
anat = anii.img;
anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
anat = anat - maxt;
anat2 = anat; anat2(p2tmp < Pthresh) = ttmp(p2tmp < Pthresh); anat2(mask == 0) = 0;
anat1 = anat; anat1(p1tmp < Pthresh) = ttmp(p1tmp < Pthresh); anat1(mask == 0) = 0;
for i = 1:size(nii.img,3)
    figure(h2); subplot(6,4,i);
    imagesc(p2tmp(:,:,i)); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbarp); caxis([0 1.1])
    
    figure(h3); subplot(6,4,i); 
    imagesc(p1tmp(:,:,i)); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbarp); caxis([0 1.1])
    
    figure(h4); subplot(6,4,i);
    imagesc(anat2(:,:,i)); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbar); caxis(CAX); 
    
    figure(h5); subplot(6,4,i);
    imagesc(anat1(:,:,i)); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbar); caxis(CAX); 
    
end

figure(h1); title('GLM')
figure(h2); title(['2 tail - P values'])
figure(h3); title(['1 tail - P values'])
figure(h4); title(['2 tail - P < ' num2str(Pthresh)])
figure(h5); title(['1 tail - P < ' num2str(Pthresh)])

if isequal(template,1)
    t = '_template';
    func_fold = [func_fold 'template' slash];
else
    t = [];
    func_fold = [func_fold 'notemplate' slash];
end
if ~exist(func_fold,'dir'); mkdir(func_fold); end
savefig(h1,[func_fold 'Design.fig'])
savefig(h2,[func_fold 'pval_2tail' t '.fig'])
savefig(h3,[func_fold 'pval_1tail' t '.fig'])
savefig(h4,[func_fold 'tval_p2tail' t '.fig'])
savefig(h5,[func_fold 'tval_p1tail' t '.fig'])

saveas(h1,[func_fold 'Design.png'])
saveas(h2,[func_fold 'pval_2tail' t '.png'])
saveas(h3,[func_fold 'pval_1tail' t '.png'])
saveas(h4,[func_fold 'tval_p2tail' t '.png'])
saveas(h5,[func_fold 'tval_p1tail' t '.png'])
