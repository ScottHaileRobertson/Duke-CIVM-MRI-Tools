% clc; close all;
% 
% load('mri'); % Matlab default
% D = double(squeeze(D));
% 
% % Show volume before filtration
% figure();
% imslice(abs(D),'Unfiltered volume');
% 
% % Bilateral filter parameters
% kernel_width = 5; % how far the filter filters each voxel
% % intensity_sigma = 0.003; % a rough measurement of intensity similarity
% intensity_sigma = 0.002; % a rough measurement of intensity similarity
% spatial_sigma = 1.25; % a rough measurement of spatial closeness
% 
% % Apply bilateral filter
% c = BF_3D(abs(b),kernel_width,intensity_sigma,spatial_sigma);
% 
% % Show filtered volume
% figure();
% imslice(abs(c),'Filtered volume');

clc; close all;

% Bilateral filter parameters
kernel_radius = 7; % Size of filter kernel - Darin says to use at least 13 to have accurate DC info
intensity_sigma = 0.4; % Intensity filtration (domain)
spatial_sigma = 1.2; % spatial_filtration (range)

% Create noisy phantom
matrixSize = 128;
nSlices = 128;
noiseFloor = 0.2;
noisyVol = phantom(matrixSize,matrixSize);
noisyVol = repmat(noisyVol,[1 1 nSlices]);
noisyVol = noisyVol + noiseFloor*randn([matrixSize,matrixSize,nSlices]);

% Mirror the input image to handle edges
noisyVol = padarray(noisyVol,kernel_radius*[1 1 1],'symmetric','both');

% Apply bilateral filter
tic;
disp('Bilateral Filtering...');
filtVol = BF_3D(noisyVol,kernel_radius,intensity_sigma,spatial_sigma);
toc

% Crop mirroring
filtVol = filtVol((kernel_radius+1):(end-kernel_radius),...
	(kernel_radius+1):(end-kernel_radius),...
	(kernel_radius+1):(end-kernel_radius));
noisyVol = noisyVol((kernel_radius+1):(end-kernel_radius),...
	(kernel_radius+1):(end-kernel_radius),...
	(kernel_radius+1):(end-kernel_radius));

% Show volume before filtration
figure();
imslice(noisyVol,'Unfiltered volume');

% Show filtered volume
figure();
imslice(filtVol,'Filtered volume');

% Show difference volume
figure();
imslice(noisyVol - filtVol,'Difference volume');