%% Load Volumes
disp('Load Ventilation volume');
vent_path = filepath('/home/scott/Public/data/')
load(vent_path);

disp('Load Proton volume');
proton_path = filepath('/home/scott/Public/data/')
load(proton_path);

%% Register volumes
iFixed = abs(gasVol);
iMoving = abs(uteVol);

% Add info about two volumes
Rfixed  = imref3d(size(iFixed),1,1,1);
Rmoving = imref3d(size(iMoving),1,1,1);

% Set up initial registration with a rigid transformation (fast, but not
% super accurate
[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.004;
optimizer.MaximumIterations=300;
disp('Starting rough rigid registration');
iRegistered = imregister(iMoving,Rmoving, iFixed,Rfixed, 'rigid', optimizer, metric);
Rmoving = imref3d(size(iRegistered),1,1,1);

% Perform a better affine transformation
optimizer.Epsilon = 1E-40;
optimizer.MaximumIterations=500;
disp('Starting fine affine registration');
iRegistered = imregister(iRegistered,Rmoving, iFixed,Rfixed, 'affine', optimizer, metric);

% Show registered
half_size = round(0.5*size(iFixed));
subplot(3,2,1);
imshowpair(squeeze(iFixed(half_size(1),:,:)),squeeze(iMoving(half_size(1),:,:)));
title('Pre');
subplot(3,2,3);
imshowpair(squeeze(iFixed(:,half_size(2),:)),squeeze(iMoving(:,half_size(2),:)));
subplot(3,2,5);
imshowpair(squeeze(iFixed(:,:,half_size(3))),squeeze(iMoving(:,:,half_size(3))));

subplot(3,2,2);
imshowpair(squeeze(iFixed(half_size(1),:,:)),squeeze(iRegistered(half_size(1),:,:)));
title('Post');
subplot(3,2,4);
imshowpair(squeeze(iFixed(:,half_size(2),:)),squeeze(iRegistered(:,half_size(2),:)));
subplot(3,2,6);
imshowpair(squeeze(iFixed(:,:,half_size(3))),squeeze(iRegistered(:,:,half_size(3))));

[pathstr,name,ext] = fileparts(proton_path);
save([pathstr name '_reg' '.mat'],'iRegistered');



