function [fix,float,coreg,other,otherc] = fUS_Coreg(fix,float,other)


float = imresize(float,size(fix));

% Type of registration error used see registration_error.m
type='sd';

% Smooth both images for faster registration
I1s=imfilter(float,fspecial('gaussian'));
I2s=imfilter(fix,fspecial('gaussian'));

% Parameter scaling of translateX translateY rotate resizeX resizeY 
% (ShearXY, ShearYX can also be included add the end of the vector, but 
%   because rotation can be seen as a type of shear transform, the scaling 
%       of the rotation or shear transform must be very small)
scale=[1 1 1 0.01 0.01];

[x]=lsqnonlin(@(x)affine_registration_image(x,scale,I1s,I2s,type),[0 0 0 100 100],[],[],optimset('Display','iter','MaxIter',100));

% Scale the translation, resize and rotation parameters to the real values
x=x.*scale;

% Make the affine transformation matrix
M=make_transformation_matrix(x(1:2),x(3),x(4:5));

coreg=affine_transform(float,M,3); % 3 stands for cubic interpolation

% Also functional images?
if ~isempty(other)
    
    for i = 1:size(other,3)
        otherc(:,:,i) = affine_transform(other(:,:,i),M,3);
        
        if size(other,3)<10
            subplot(size(other,3),3,1+3*(i-1)); imagesc(other(:,:,i));
            title('orig'); set(gca,'xticklabel',[],'yticklabel',[]);
            subplot(size(other,3),3,2+3*(i-1)); imagesc(otherc(:,:,i));
            title('coreg'); set(gca,'xticklabel',[],'yticklabel',[])
            subplot(size(other,3),3,3+3*(i-1)); imagesc(other(:,:,i) - otherc(:,:,i));
            title('orig - coreg'); set(gca,'xticklabel',[],'yticklabel',[])
        end
    end
end

% Show the registration results
figure(76)
subplot(2,3,1), imagesc(float); title('float'); set(gca,'xticklabel',[],'yticklabel',[])
subplot(2,3,2), imagesc(fix); title('fix'); set(gca,'xticklabel',[],'yticklabel',[])
subplot(2,3,3), imagesc(coreg); title('float coreg'); set(gca,'xticklabel',[],'yticklabel',[])
subplot(2,3,4), imagesc(fix - float); title('fix - float'); set(gca,'xticklabel',[],'yticklabel',[])
subplot(2,3,5), imagesc(fix - coreg); title('fix - coreg'); set(gca,'xticklabel',[],'yticklabel',[])
subplot(2,3,6), imagesc(float - coreg); title('float - coreg'); set(gca,'xticklabel',[],'yticklabel',[])
pause(5)
