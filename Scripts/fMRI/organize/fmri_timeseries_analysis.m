function fmri_timeseries_analysis(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
descr = param.descr;

% Load group data
for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_grp_data.mat'];
    count_file = [stim_storage stim{i_stim} '_active_voxel_count.mat'];
    save_file = [stim_storage stim{i_stim} '_time_series_info.mat'];
    
    if exist(stim_file,'file') && exist(count_file,'file')
        
        load(stim_file)
        load(count_file)
        
        % Time series analysis for each ROI
        start = [30 70 110 150 190]; % Stimulus start times
        [b,a] = butter(4,[0.006 0.1]/((1/1.5)/2)); % time-domain filter
        for i_mouse = 1:size(base_fold,1)
            
            S = size(grp_data);
            datatmp = reshape(grp_data(:,:,:,i_mouse),S(1)*S(2),S(3));
            
            for i_pix = 1:size(datatmp,1)
                
                ts{i_pix,i_mouse} = datatmp(i_pix,:); % sum time series in roi
                % baseline normalize
                meanbase = mean(ts{i_pix,i_mouse}(1:20)); % baseline of time series
                tsnorm{i_pix,i_mouse} = (ts{i_pix,i_mouse} - meanbase)/meanbase*100; % baseline normalized roi time series
                
                % filter normalized time series
                tmp = [zeros(1,100) tsnorm{i_pix,i_mouse} zeros(1,100)];
                tmp = filtfilt(b,a,tmp);
                tsnorm{i_pix,i_mouse} = tmp(101:end-100); % filtered normalized time series 
                
                % also filter un-normalized time series for z score compute
%                 tmp = [zeros(1,100) ts{i_pix,i_mouse} zeros(1,100)];
                ts{i_pix,i_mouse} = filtfilt(b,a,ts{i_pix,i_mouse});
%                 ts{i_pix,i_mouse} = tmp(101:end-100);
                
                ts_normpeaks = zeros(size(start,2)-1,size(-10:28,2));
                for k = 1:size(start,2)
                    ts_normpeaks(k,:) = tsnorm{i_pix,i_mouse}(start(k)-10:start(k)+28);
                end
                ts_normpeaks = ts_normpeaks'-repmat(mean(ts_normpeaks',1),[39 1]);
                ts_normpeaks = ts_normpeaks - repmat(mean(ts_normpeaks(1:10,:),1),[39,1]);
                
                % z-score conversion
                base = 1:20;
                tmp = (ts{i_pix,i_mouse} - mean(ts{i_pix,i_mouse}(base)))/...
                    std(ts{i_pix,i_mouse}(base));
                for k = 1:size(start,2)
                    z(k) = mean(tmp(start(k):start(k)+8));
                end
%                 tsnormz = (tmp - mean(tmp(base)))/std(tmp(base));
                
                tsnormave{i_pix,i_mouse} = mean(ts_normpeaks',1);
                tsnormavepeak{i_pix,i_mouse} = max(mean(ts_normpeaks'));
                tsnormaveAUC{i_pix,i_mouse} = sum(mean(ts_normpeaks'));
                tsnormZ{i_pix,i_mouse} = tmp;
                tsnormaveZ{i_pix,i_mouse} = mean(z);
                
            end
            

            % Create noise roi
            if ~exist('noisepts','var')
                tmp = zeros(35,80);
                allpts = vertcat(pts{:,1});
                tmp(allpts) = 1;
                f = figure; imagesc(tmp)
                [xnoise, ynoise] = getpts(f);
                t = linspace(0, 2*pi); r = 3;
                xnoise = r * cos(t) + xnoise; ynoise = r * sin(t) + ynoise;
                c2 = patch(xnoise, ynoise, zeros(1,size(xnoise,2)),...
                    'linewidth',1.5,'facecolor','none','edgecolor','g');
                tmpnoise = unique([round(xnoise)' round(ynoise)'],'rows');
                tmpnoise((sum(tmpnoise<=0,2)>0),:) = [];
                tmp(tmpnoise(:,2),tmpnoise(:,1)) = 2;
                figure(f); cla; imagesc(tmp)
                noisepts = find(tmp == 2);
                close(f)
            end
            
            
            for i_roi = 1:size(pts,1) + 1
                
                if i_roi == size(pts,1) + 1
                    PTS = noisepts;
                else
                    PTS = pts{i_roi,i_mouse};
                end
                
                ts_roi{i_roi,i_mouse} = sum(datatmp(PTS,:),1); % sum time series in roi
                meanbase = mean(ts_roi{i_roi,i_mouse}(1:20)); % baseline of time series
                tsnorm_roi{i_roi,i_mouse} = (ts_roi{i_roi,i_mouse} - meanbase)/meanbase*100; % baseline normalized roi time series

                p = polyfit(1:100,tsnorm_roi{i_roi,i_mouse}(1:100),1);
                if p(1) > 0
                    tmp1 = [-tsnorm_roi{i_roi,i_mouse}(1:100)+tsnorm_roi{i_roi,i_mouse}(100)];
                else
                    tmp1 = [tsnorm_roi{i_roi,i_mouse}(1:100)-tsnorm_roi{i_roi,i_mouse}(100)];
                end
                
                p = polyfit(1:100,tsnorm_roi{i_roi,i_mouse}(end-99:end),1);
                if p(1) > 0
                    tmp2 = [tsnorm_roi{i_roi,i_mouse}(end-99) - tsnorm_roi{i_roi,i_mouse}(end-99:end)];
                else
                    tmp2 = [tsnorm_roi{i_roi,i_mouse}(end-99) + tsnorm_roi{i_roi,i_mouse}(end-99:end)];
                end
                tmp = [tmp1 tsnorm_roi{i_roi,i_mouse} tmp2];

                %                 tmp = [max(tsnorm_roi{i_roi,i_mouse})*ones(1,100) tsnorm_roi{i_roi,i_mouse} min(tsnorm_roi{i_roi,i_mouse})*ones(1,100)];
                tmp = filtfilt(b,a,tmp);
                tsnorm_roi{i_roi,i_mouse} = tmp(101:end-100);
                
%                 p = polyfit(1:100,ts_roi{i_roi,i_mouse}(1:100),1);
%                 if p(1) > 0
%                     tmp1 = [-ts_roi{i_roi,i_mouse}(1:100)+ts_roi{i_roi,i_mouse}(100)];
%                 else
%                     tmp1 = [ts_roi{i_roi,i_mouse}(1:100)-ts_roi{i_roi,i_mouse}(100)];
%                 end
%                 
%                 p = polyfit(1:100,ts_roi{i_roi,i_mouse}(end-99:end),1);
%                 if p(1) > 0
%                     tmp2 = [ts_roi{i_roi,i_mouse}(end-99) - ts_roi{i_roi,i_mouse}(end-99:end)];
%                 else
%                     tmp2 = [ts_roi{i_roi,i_mouse}(end-99) + ts_roi{i_roi,i_mouse}(end-99:end)];
%                 end
%                 tmp = [tmp1 ts_roi{i_roi,i_mouse} tmp2];
%                 
%                 % also filter un-normalized time series for z score compute
%                 tmp = [zeros(1,100) ts_roi{i_roi,i_mouse} - mean(ts_roi{i_roi,i_mouse}) zeros(1,100)];
                ts_roi{i_roi,i_mouse} = filtfilt(b,a,ts_roi{i_roi,i_mouse});
%                 ts_roi{i_roi,i_mouse} = tmp(101:end-100);
                
                % tsnorm_roi{i_roi,i_mouse} = filtfilt(b,a,tsnorm_roi{i_roi,i_mouse}); % filtered normalized time series
                
                ts_normpeaks = zeros(size(start,2)-1,size(-10:28,2));
                for k = 1:size(start,2)
                    ts_normpeaks(k,:) = tsnorm_roi{i_roi,i_mouse}(start(k)-10:start(k)+28);
                end
                ts_normpeaks = ts_normpeaks'-repmat(mean(ts_normpeaks',1),[39 1]);
                ts_normpeaks = ts_normpeaks - repmat(mean(ts_normpeaks(1:10,:),1),[39,1]);
                
                % z-score conversion
                base = 1:20;
                ts_roi{i_roi,i_mouse} = detrend(ts_roi{i_roi,i_mouse});
                tmp = (ts_roi{i_roi,i_mouse} - mean(ts_roi{i_roi,i_mouse}(base)))/...
                    std(ts_roi{i_roi,i_mouse}(base));
                for k = 1:size(start,2)
                    z(k) = mean(tmp(start(k):start(k)+8));
                end
%                 tsnormz = (tmp - mean(tmp(base)))/std(tmp(base));
                
                if i_roi == size(pts,1) + 1
                    tsnormave_noise{i_mouse} = mean(ts_normpeaks',1);
                    tsnormavepeak_noise{i_mouse} = max(mean(ts_normpeaks'));
                    tsnormaveAUC_noise{i_mouse} = sum(mean(ts_normpeaks'));
                    tsnormZ_noise{i_roi,i_mouse} = tmp;
                    tsnormaveZ_noise{i_roi,i_mouse} = mean(z);
                else
                    tsnormave_roi{i_roi,i_mouse} = mean(ts_normpeaks',1);
                    tsnormavepeak_roi{i_roi,i_mouse} = max(mean(ts_normpeaks'));
                    tsnormaveAUC_roi{i_roi,i_mouse} = sum(mean(ts_normpeaks'));
                    tsnormZ_roi{i_roi,i_mouse} = tmp;
                    tsnormaveZ_roi{i_roi,i_mouse} = mean(z);
                end
                
            end
            
        end
        
    end
    
    fprintf('\nSaving: %s\n', save_file);
    save(save_file,'tsnorm','tsnormave','tsnormavepeak','tsnormaveAUC',...
        'tsnormZ','tsnormaveZ','tsnormZ_roi','tsnormaveZ_roi','tsnormZ_noise','tsnormaveZ_noise',...
        'tsnorm_roi','tsnormave_roi','tsnormavepeak_roi','tsnormaveAUC_roi',...
        'tsnormave_noise','tsnormavepeak_noise','tsnormaveAUC_noise');
    cell2mat(tsnormaveZ_noise(end,:))
end         
                
x = 1;
                
                
                
                