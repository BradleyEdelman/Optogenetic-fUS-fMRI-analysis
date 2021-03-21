function fus_organize_data(storage,base_fold,slash,param)

stim = {'0_1' '0_5' '1_0'};
stim1 = {'0.1','0.5','1.0'};
Dummy = param.Dummy;
descr = param.descr;

% Load individual animal data and construct group-level design matrix
for i_stim = 1:size(stim,2)
    
    clear ind_data grp_data glmstim glmfilt grp_glm
    grp_data = []; D1 = []; D2 = []; D3 = [];
    for i_mouse = 1:size(base_fold,1)

        % Specify reconstructed functional scans from session
        recon_fold = [storage base_fold{i_mouse}(1:8) slash base_fold{i_mouse}(10:end) slash];
        snr_fold = [recon_fold stim1{i_stim} slash];
        func_file = [snr_fold 'template' slash stim1{i_stim} '_Func.mat'];

        if exist(func_file,'file')
            
            % Load and store individual data file
            load(func_file)
            ind_data{i_mouse} = imresize3(FUNC.Data.pdi(:,:,Dummy+1:end), [110 80 230]);
            
            % Extract components of individual GLM to creat group GLM
            glmstim{i_mouse} = FUNC.GLM.X(:,1:5);
            glmfilt{i_mouse} = FUNC.GLM.X(:,6:12);
            D1 = blkdiag(D1,glmstim{i_mouse}(:,1:4));
            D2 = blkdiag(D2,glmstim{i_mouse}(:,5));
            D3 = blkdiag(D3,glmfilt{i_mouse});
            
            % Compile group data
            grp_data = cat(4,grp_data,...
                imresize3(FUNC.Data.pdi(:,:,Dummy+1:end), [110 80 230]));

        end

    end
    grp_glm = [D1 D2 D3];
    
    % Save individual and group data
    stim_storage = [storage stim{i_stim} slash descr slash];
    if ~exist(stim_storage,'dir'); mkdir(stim_storage); end
    save_file = [stim_storage stim{i_stim} '_grp_data.mat'];
    fprintf('\nSaving: %s\n', save_file);
    save(save_file,'grp_data','glmstim','glmfilt','grp_glm')
    
end

    
    
    
    
    
    