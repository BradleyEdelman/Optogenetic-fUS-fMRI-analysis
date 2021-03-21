fixed = mat2gray(pdi(:,:,1));
moving = mat2gray(pdi(:,:,2));
%%

fixed = pdi(:,:,1); fixed = 10*log10(fixed./max(max(max(pdi)))); fixed = mat2gray(fixed);
for i = 2:100%size(pdi,3)
    i
    moving = pdi(:,:,i);
    moving = 10*log10(moving./max(max(max(pdi))));
    moving = mat2gray(moving);
    
    tform = imregtform(moving,fixed,'rigid',optimizer,metric);
    movingReg(:,:,i) = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
    
    dx(i) = tform.T(3,1);
    dy(i) = tform.T(3,2);
    dtheta(i) = tform.T(1,2);
    
end

figure;
subplot(3,1,1); plot(dx)
subplot(3,1,2); plot(dy)
subplot(3,1,3); plot(dtheta)

%%
figure, imshowpair(fixed,moving,'montage');
tformEstimate = imregcorr(moving,fixed,'rigid');

movingReg = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));


figure, imshowpair(fixed,movingReg,'montage');
%%








fixed  = imread('cameraman.tif');
theta = 0;
S = 1.0;
tx = 5;
ty = 5;
tformFixedToMoving = affine2d([S.*cosd(theta) -S.*sind(theta) 0; S.*sind(theta) S.*cosd(theta) 0; tx ty 1]);
moving = imwarp(fixed,tformFixedToMoving,'OutputView',imref2d(size(fixed)));
% Add a bit of noise uniform noise to moving to make the problem harder.
% moving = moving + uint8(10*rand(size(moving)));
tformMovingToFixedEstimate = imregcorr(moving,fixed);
figure, imshowpair(fixed,moving,'montage');
% Apply estimated geometric transform to moving. Specify 'OutputView' to
% get registered moving image that is the same size as the fixed image.
Rfixed = imref2d(size(fixed));
movingReg = imwarp(moving,tformMovingToFixedEstimate,'OutputView',Rfixed);
figure, imshowpair(fixed,movingReg,'montage');
tformFixedToMovingEst = invert(tformMovingToFixedEstimate);
tformFixedToMovingEst.T
angleErrorInDegrees = acosd(tformFixedToMovingEst.T(1,1)) - theta
translationError = [tx ty] - tformFixedToMovingEst.T(3,1:2)




