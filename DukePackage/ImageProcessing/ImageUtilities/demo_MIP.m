clc; clear all; close all;
load('thresh_vol.mat');

b = c;
clear c;
angles = 1:1:360;
numAngles = length(angles);
dims = size(b);
rotVol = zeros(dims);
mipImages = zeros(dims(1), dims(2), numAngles);
for i = 1:numAngles
    ang = angles(i)
    
    nSlices = size(b,1);
    
    %Create rotated volume
    for j = 1:nSlices
        rotVol(j,:,:) = imrotate(squeeze(b(j,:,:)), ang, 'bicubic', 'crop');
%         rotVol(j,:,:) = imrotate(squeeze(b(j,:,:)), ang, 'bilinear', 'crop');
    end
    
    %Create MIP image
    mipImages(:,:,i) = max(rotVol,[],3);
%     mipImages(:,:,i) = sum(rotVol,3);
end

% Prepare the new file.
vidObj = VideoWriter('Recon_MIP_bicubic');
vidObj.FrameRate = 30;
open(vidObj);
fig = figure();

% Create movie
minmax = [min(mipImages(:)) max(mipImages(:))];
for i=1:numAngles
    imagesc(mipImages(:,:,i),minmax);
    colormap(gray);
    axis image;
    
    writeVideo(vidObj,getframe(fig));
end

% Close the file.
close(vidObj);

% Show the mip
imslice(mipImages);