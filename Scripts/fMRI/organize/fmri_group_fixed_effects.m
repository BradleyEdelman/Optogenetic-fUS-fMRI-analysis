function fmri_group_fixed_effects(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
descr = param.descr;

% Define t-ditributions
tdist2T = @(t,DF) (1-betainc(DF/(DF+t^2),DF/2,0.5)); % 2-tailed

% Load group data
for i_stim = 1:size(stim,2)
    
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    stim_file = [stim_storage stim{i_stim} '_grp_data.mat'];
    save_file = [stim_storage stim{i_stim} '_grp_fixed_effects.mat'];
    
    if exist(stim_file,'file')
        
        load(stim_file)
        
        DF = size(grp_glm,1) - 2;
        grp_D = grp_glm;
        grp_c = [repmat([1 0 0 0],[1 size(base_fold,1)]),...
            zeros(1,size(base_fold,1)),zeros(1,7*size(base_fold,1))]';
        
        for i = 1:size(grp_data,1)
            for j = 1:size(grp_data,2)
                
                Y = squeeze(grp_data(i,j,:));
                grp_beta_hat = inv(grp_D'*grp_D)*grp_D'*Y;
                    
                grp_beta(i,j,:) = reshape(grp_beta_hat,[1 1 size(grp_c,1)]);

                grp_Var_e = (Y-grp_D*grp_beta_hat)'*(Y-grp_D*grp_beta_hat)/DF;

                %Hypothesis testing; Compute the t statistic
                grp_t_stat(i,j)=grp_c'*grp_beta_hat/sqrt(grp_Var_e*grp_c'*inv(grp_D'*grp_D)*grp_c);

                if isnan(grp_t_stat(i,j))
                    grp_p2tail(i,j) = 1;
                else
                    grp_p2tail(i,j) = 1 - tdist2T(grp_t_stat(i,j),DF);
                end
                
            end
        end
        
        [h, crit_p, adj_ci_cvrg, grp_c_p2tail]=fdr_bh(grp_p2tail,param.Pthresh,'pdep','no');
        
        fprintf('\nSaving: %s\n', save_file);
        save(save_file,'grp_beta','grp_t_stat','grp_p2tail','grp_c_p2tail','grp_c','grp_D');
        
    end
end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    