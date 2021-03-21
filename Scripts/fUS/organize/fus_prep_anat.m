function [anat,Cbar,CAX] = fus_prep_anat(storage,slash,param,plotmin,plotmax,MAX)

% Create colorbar for activation maps
% plotmin = -100; plotmax = 0; maxt = param.maxt; %Pthresh = param.Pthresh;
Cbar = [gray(abs(plotmax - plotmin)*100); jet(2*MAX*100)]; CAX = [plotmin - MAX MAX];
        

% Load anatomical image
if isfield(param,'fusroi')
    atmp = [storage '20191119' slash '4346075_N_D2' slash 'Anat.mat'];
    anii = load(atmp); anii = anii.ref;
    mask = ones(size(anii));
elseif isfield(param,'fmriroi')
    anii = load_nii(param.atmp); anii = anii.img(:,:,10);
    mnii = load_nii(param.mtmp); mask = mnii.img(:,:,10);
end

 % rescale dynamic range of anatomical for visuslization
anat = anii;
anat = (anat - min(anat(:)))/(max(anat(:)) - min(anat(:)))*(plotmax - plotmin) + plotmin;
anat = anat - MAX;