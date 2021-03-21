function fmri_compute_snr(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
descr = param.descr;

% Load group data
for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_grp_data.mat'];
    count_file = [stim_storage stim{i_stim} '_active_voxel_count.mat'];
    
    if exist(stim_file,'file') && exist(count_file,'file')
        
        load(stim_file)
        load(count_file)
        
        Start = [30 70 110 150 190];
        End = [38 78 118 158 198];
        for i_mouse = 1:size(base_fold,1)
            
            Data = grp_data(:,:,:,i_mouse);
%             Data = Data(:,:,1:Start(1));
%             Data(:,:,[Start' End']) = [];
            figure(1); subplot(1,size(base_fold,1),i_mouse); imagesc(Data(:,:,1));
            if isequal(i_stim,4)
                f=figure; imagesc(Data(:,:,1));  axis tight; colormap gray
                [xsig{i_mouse}, ysig{i_mouse}] = getpts(f);
                [xnoise{i_mouse}, ynoise{i_mouse}] = getpts(f);
             
                t = linspace(0, 2*pi); r = .5;

                xsig{i_mouse} = r * cos(t) + xsig{i_mouse}; ysig{i_mouse} = r * sin(t) + ysig{i_mouse};
                c1(i_mouse) = patch(xsig{i_mouse}, ysig{i_mouse}, zeros(1,size(xsig{i_mouse},2)),...
                    'linewidth',1.5,'facecolor','none','edgecolor','r');
                xnoise{i_mouse} = r * cos(t) + xnoise{i_mouse}; ynoise{i_mouse} = r * sin(t) + ynoise{i_mouse};
                c2(i_mouse) = patch(xnoise{i_mouse}, ynoise{i_mouse}, zeros(1,size(xnoise{i_mouse},2)),...
                    'linewidth',1.5,'facecolor','none','edgecolor','g');
            end
            
            RM1 = 14;
            pts2 = vertcat(pts{[13:14],i_mouse});
            pts2 = pts{RM1,i_mouse};
            sz = [size(Data,1) size(Data,2)];
            [row,col] = ind2sub(sz,pts2);
            sig = Data(row,col,:);
            noise = Data(5:9,65:69);
            
            snr(i_mouse,i_stim) = mean2(sig(:,:,1))/std2(noise(:,:,1));
            tsnr(i_mouse,i_stim) = mean(sum(sum(sig,1),2))/std(sum(sum(sig,1),2));
            
%             tmpsig = unique([round(xsig{i_mouse})' round(ysig{i_mouse})'],'rows');
%             tmpnoise = unique([round(xnoise{i_mouse})' round(ynoise{i_mouse})'],'rows');
%             tmpsig((sum(tmpsig<=0,2)>0),:) = [];
%             tmpnoise((sum(tmpnoise<=0,2)>0),:) = [];
%             
%             % snr
%             sigInt = Data(tmpsig(:,2),tmpsig(:,1),1);
%             noiseInt = Data(tmpnoise(:,2),tmpnoise(:,1),1);
%             snr(i_mouse,i_stim) = mean2(sigInt)/std2(noiseInt);
% 
%             % tsnr
%             sigInt = squeeze(sum(sum(Data(tmpsig(:,2),tmpsig(:,1),:),1),2));
%             noiseInt = squeeze(sum(sum(Data(tmpnoise(:,2),tmpnoise(:,1),:),1),2));
%             tsnr(i_mouse,i_stim) = mean(sigInt)/std(sigInt);
%             
        end
    end
    
end
save_file = [storage descr '_SNR.mat'];
fprintf('\nSaving: %s\n', save_file);
save(save_file,'snr','tsnr');

h1 = figure; hold on
snr = mean(snr,2); tsnr = mean(tsnr,2);
Msnr = mean(snr); Mtsnr = mean(tsnr);
Ssnr = std(snr)/sqrt(size(snr,1)); Stsnr = std(tsnr)/sqrt(size(snr,1));

bar([Msnr Mtsnr])
errorbar([Msnr Mtsnr],[Ssnr Stsnr]);
set(gca,'xtick',[1 2],'xticklabel',{'SNR' 'tSNR'},'ylim',[0 1000])
title([descr 'SNR'])
saveas(h1,[storage descr '_SNR.svg']);
saveas(h1,[storage descr '_SNR.fig']);
            
            
            
            
            
            
            
            
            
            