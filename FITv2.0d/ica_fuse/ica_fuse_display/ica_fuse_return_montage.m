function [im] = ica_fuse_return_montage(images, DIM)
% returns the montage of the image
% 
% Input:
% 1. images - 
% 
% Ouput:
% 1. im - montage


if length(DIM) == 2
    DIM(3) = 1;
end

% --get montage of image
images = squeeze(images);
a1b = zeros([DIM(2), DIM(1), 1, DIM(3)]);

temp = reshape(images, [DIM(1), DIM(2), 1, DIM(3)]);

clear images;
for k=1:DIM(3)
    a1b(:,:,1,end-k+1) = (temp(:,end:-1:1,1,k)');
end
clear temp;
temp=reshape(a1b(:,:,:,:), [DIM(2), DIM(1), 1, DIM(3)]);
[im] = ica_fuse_get_montage(temp, -999);