
post = load_nii('D:\Data_Processed\fMRI\20191124\BEd_postCW_4364122_N_1_4\0.1\0.1_EPI.nii');
postm = load_nii('D:\Data_Processed\fMRI\20191124\BEd_postCW_4364122_N_1_4\13a\ave_anatomy_mask.nii');

m = postm.img(:,:,9);
m = imresize(m,[35 80]);
m = circshift(m,[-1 2]);
e = post.img(:,:,7,1);
e(m == 0) =0;
figure(2); clf
subplot(1,2,1); imagesc(e);
colormap gray

pre = load_nii('D:\Data_Processed\fMRI\20191111\BEd_preCW_4364122_N_D2_1_53\0.1\0.1_EPI.nii');
prem = load_nii('D:\Data_Processed\fMRI\20191111\BEd_preCW_4364122_N_D2_1_53\18a\18_anatomy_mask.nii');

m = prem.img(:,:,9);
m = imresize(fliplr(rot90(m)),[35 80]);
m = imresize(m,0.9);
m = padarray(m,[35 - size(m,1) 80 - size(m,2)],0);
m = m(1:35,1:80);
m = imresize(m,2);
m = circshift(m,[-1 -10]);
m = imresize(m,.5);
e = pre4.img(:,:,9,1);
e(m == 0) =0;
subplot(1,2,2); imagesc(e);
colormap gray
