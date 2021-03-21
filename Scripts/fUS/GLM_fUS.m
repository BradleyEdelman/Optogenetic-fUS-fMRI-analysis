function [Block,GLM] = GLM_fUS(func_fold,anat_fold,fpdi,param)

if isunix; slash = '/'; elseif ispc; slash = '\'; end

dummy = param.Dummy;
order = param.order;
template = param.template;
templateidx = param.templateidx;

if ~isempty(anat_fold)
    
    if isnan(templateidx)
        ave_anat = [anat_fold 'Ave_Anat.mat'];
    else
        ave_anat = [anat_fold 'Total_Anat_orig.mat'];
    end

    if exist(ave_anat,'file') && isequal(template,1)
        
        tmp = load(ave_anat);
        if isnan(templateidx)
            apdi = tmp.aveimg;
        else
            apdi = tmp.totimg_orig(:,:,templateidx);
        end
        
        ave_anat_mask = [anat_fold 'Ave_Anat_Mask.mat'];
        if exist(ave_anat_mask,'file')
            tmp = load(ave_anat_mask);
            mask = tmp;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            mask = ones(size(apdi));
        end
        
    else
        
        anat = [anat_fold 'Anat.mat'];
        tmp = load(anat);
        
        if isfield(tmp,'source')
            apdi = tmp.source;
        elseif isfield(tmp,'ref')
            apdi = tmp.ref;
        end
        
        anat_mask = [anat_fold 'Anat_Mask.mat'];
        if exist(anat_mask,'file')
            tmp = load(anat_mask);
            mask = tmp;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            mask = ones(size(apdi));
        end
    end
    
else

    apdi = fpdi(:,:,1);
    mask = ones(size(apdi));
        
end


Block = [];
if isunix
    load('/media/bradley/Seagate Backup Plus Drive/Preprocessing/4th_order_gam_Design_fUS.mat')
elseif ispc
    load('D:\Preprocessing\4th_order_gam_Design_fUS.mat')
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
fpdi = fpdi(:,:,dummy+1:end);
nTRs = size(fpdi,3);

% Define t-ditributions
tdist2T = @(t,DF) (1-betainc(DF/(DF+t^2),DF/2,0.5)); % 2-tailed
tdist1T = @(t,DF) 1-(1-tdist2T(t,DF))/2; % 1-tailed

c = [1;zeros(size(X,2)-1,1)];
beta = zeros(size(fpdi,1),size(fpdi,2),size(c,1));
t_stat = zeros(size(fpdi,1),size(fpdi,2));
p2tail = ones(size(fpdi,1),size(fpdi,2));
p1tail = ones(size(fpdi,1),size(fpdi,2));
for i=1:size(fpdi,1)
    for j=1:size(fpdi,2)

        Y = squeeze(fpdi(i,j,:))';
        
        DF = nTRs - 2;
        
        beta_hat=inv(X'*X)*X'*Y';
        beta(i,j,:) = reshape(beta_hat,1,1,size(c,1));
        Var_e=(Y'-X*beta_hat)'*(Y'-X*beta_hat)/DF;
        
        %Hypothesis testing; Compute the t statistic
        t_stat(i,j)=c'*beta_hat/sqrt(Var_e*c'*inv(X'*X)*c);
        
        if isnan(t_stat(i,j))
            p2tail(i,j) = 1; 
            p1tail(i,j) = 1;
        else
            p2tail(i,j) = 1- tdist2T(t_stat(i,j),DF);
            p1tail(i,j) = 1- tdist1T(t_stat(i,j),DF);
        end

    end
end

GLM.X =X;
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

h2 = figure(11);
for i = 1:size(c,1)
    subplot(1,size(c,1),i); imagesc(beta(:,:,i));
    set(gca,'yticklabel',[],'xticklabel',[])
    AX(i,:) = caxis;
end
for i = 1:size(c,1); subplot(1,size(c,1),i); caxis([min(AX(:,1)) max(AX(:,2))]); end
colormap jet

Pthresh = 0.05; maxt = 50;

plotmin = -100; plotmax = 0;
Cbar = [gray(abs(plotmax - plotmin)*100); fireice(2*maxt*100)]; CAX = [plotmin - maxt maxt];
Cbarp = [jet(round(1.1*100));[0 0 0]];
h3 = figure(12); h4 = figure(13);


% Resize time-series param to anat sampling
p2tmp = imresize(p2tail,[size(apdi,1) size(apdi,2)],'bilinear');
p2tmp(mask == 0) = 1.1;
p1tmp = imresize(p1tail,[size(apdi,1) size(apdi,2)],'bilinear');
p1tmp(mask == 0) = 1.1;
ttmp = imresize(t_stat,[size(apdi,1) size(apdi,2)]);
ttmp(mask == 0) = 0;

% Apply correction
[h, crit_p, adj_ci_cvrg, c_p2tmp]=fdr_bh(p2tmp,.05,'pdep','no');
[h, crit_p, adj_ci_cvrg, c_p1tmp]=fdr_bh(p1tmp,.05,'pdep','no');

% rescale dynamic range of anat for vis
anat = apdi; figure(50); imagesc(apdi)
anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
anat = anat - maxt;
anat2 = anat; anat2(c_p2tmp < Pthresh) = ttmp(c_p2tmp < Pthresh); %anat2(mask == 0) = 0;
anat1 = anat; anat1(c_p1tmp < Pthresh) = ttmp(c_p1tmp < Pthresh); %anat1(mask == 0) = 0;

figure(h3);
subplot(1,2,1); imagesc(c_p2tmp); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbarp); caxis([0 1.1])
subplot(1,2,2); imagesc(c_p1tmp); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbarp); caxis([0 1.1])
figure(h4)
subplot(1,2,1); imagesc(anat2); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbar); caxis(CAX)
subplot(1,2,2); imagesc(anat1); set(gca,'xticklabel',[],'yticklabel',[]); colormap(Cbar); caxis(CAX)

if isequal(template,1)
    t = '_template';
    func_fold = [func_fold 'template' slash];
else
    t = [];
    func_fold = [func_fold 'notemplate' slash];
end
if ~exist(func_fold,'dir'); mkdir(func_fold); end
savefig(h1,[func_fold 'Design.fig'])
savefig(h2,[func_fold 'beta' t '.fig'])
savefig(h3,[func_fold 'p_val' t '.fig'])
savefig(h4,[func_fold 't_val' t '.fig'])

saveas(h1,[func_fold 'Design.png'])
saveas(h2,[func_fold 'beta' t '.png'])
saveas(h3,[func_fold 'p_val' t '.png'])
saveas(h4,[func_fold 't_val' t '.png'])
