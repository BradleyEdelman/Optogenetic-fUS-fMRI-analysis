function fus_individual_fixed_effects(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
descr = param.descr;

% Define t-ditributions
tdist2T = @(t,DF) (1-betainc(DF/(DF+t^2),DF/2,0.5)); % 2-tailed
% tdist1T = @(t,DF) 1-(1-tdist2T(t,DF))/2; % 1-tailed

% Load group data
for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_grp_data.mat'];
    save_file = [stim_storage stim{i_stim} '_ind_fixed_effects.mat'];
    
    if exist(stim_file,'file')
        
        load(stim_file)
        
        for i_mouse = 1:size(base_fold,1)
            
            DF = size(grp_data,3) - 2;
            D = [glmstim{i_mouse} glmfilt{i_mouse}];
            c = [1 0 0 0 0 zeros(1,7)]';
            
            for i = 1:size(grp_data,1)
                for j = 1:size(grp_data,2)
                    
                    Y = squeeze(grp_data(i,j,:,i_mouse));
                    beta_hat = inv(D'*D)*D'*Y;
                    
                    beta(i,j,:,i_mouse) = reshape(beta_hat,[1 1 size(c,1)]);
                    
                    Var_e = (Y-D*beta_hat)'*(Y-D*beta_hat)/DF;
                    
                    %Hypothesis testing; Compute the t statistic
                    t_stat(i,j,i_mouse)=c'*beta_hat/sqrt(Var_e*c'*inv(D'*D)*c);

                    if isnan(t_stat(i,j,i_mouse))
                        p2tail(i,j,i_mouse) = 1;
                    else
                        p2tail(i,j,i_mouse) = 1 - tdist2T(t_stat(i,j,i_mouse),DF);
                    end
                    
                end
            end
        
            contrast(:,:,i_mouse) = sum((repmat(reshape(c,1,[],size(c,1)),...
                size(grp_data,1),size(grp_data,2))).*beta(:,:,:,i_mouse),3);
            
        end
        
        fprintf('\nSaving: %s\n', save_file);
        save(save_file,'beta','t_stat','p2tail','c','D','contrast');
        
    end
end



















