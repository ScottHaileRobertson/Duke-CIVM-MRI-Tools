function demo_dixon_recon_BHUTE(varargin)

if(nargin < 1 | ~exist(varargin{1}))
    disp('Select BHUTE pfile');
    bhute_pfile = filepath('/home/scott/');
else
    bhute_pfile = varargin{1};
end

if(nargin < 2)
    disp('Select dixon pfile');
    dixon_pfile = filepath('/home/scott/');
else
    dixon_pfile = varargin{2};
end

% Required parameters
% % new 128^3
% output_image_size = 64*[1 1 1];
% kernel.sharpness = 0.35;

% old 64^3
output_image_size = 64*[1 1 1];
kernel.sharpness = 0.3;

overgrid_factor = 3;
kernel.extent = 9*kernel.sharpness;
% kernel.extent = 1.5;
verbose = 1;
nPipeIter = 25;

% Human Ventilation Parameters
pfileOverride = GE.Pfile.Pfile();

pfileOverride.rdb.rdb_hdr_user1  = 0.256; % pw_gxwa
pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
pfileOverride.rdb.rdb_hdr_user44 = 1.024; % pw_gxw/1000
pfileOverride.rdb.rdb_hdr_user22 = 0.103; %toff
pfileOverride.rdb.rdb_hdr_user23 = 101; % primeplus
% pfileOverride.rdb.rdb_hdr_user23 = 137.508; % primeplus
pfileOverride.rdb.rdb_hdr_user32 = 1;

% Gradient delays
delays.x_delay = 0.000;
delays.y_delay = 0.00;
delays.z_delay = 0.000;

%Optional parameters
revision_override = [];  %Optional override if it can't be automatically read from the pfile

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(bhute_pfile);

% Convert from Pfile format
pfile = convertLegacyPfile(pfile);

% Override header values (optional)
displayPfileHeaderInfo(pfile);
pfile = overridePfile(pfile, pfileOverride);
if(verbose)
    % Display key header info
    displayPfileHeaderInfo(pfile);
end

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% Calculate trajectory for a single radial ray
radialDistance = MRI.Trajectories.Centric.Radial.calcRadialRay(pfile, delays, output_image_size);

% 	% Only keep data during gradients on
[radialDistance, pfile] = MRI.DataProcessing.removeNonReadoutSamples(radialDistance, pfile);

% Distribute rays onto 3d sphere
traj = MRI.Trajectories.Centric.Distribute.calculate3dTrajectories(radialDistance, pfile);

% Undo loopfactor
[traj, pfile] = MRI.DataProcessing.undoloopfactor(traj, pfile);

% Pick a reasonable DC sample
% dc_sample_idx = MRI.DataProcessing.calcDCsample(max(radialDistance,[],2));

% % Calculate flip angle
% MRI.DataProcessing.calcFlipAngle(pfile, dc_sample_idx);

% Calculate SNR weights
% weights = MRI.DataProcessing.calcSNRWeights(pfile, dc_sample_idx, weight_type);

% Calculate Maximum volume size for Nyquist
MRI.DataProcessing.calculateNyquistMatrixSize(radialDistance, pfile);

% Vectorize data and traj for recon
[traj, pfile] = MRI.DataProcessing.vectorizeDataAndTraj(traj, pfile);

% Enforce Nyquist limits
[traj, pfile] = MRI.DataProcessing.enforceNyquistBounds(traj, pfile);

%% Reconstruct data
% Choose kernel
kernelObj = Recon.SysModel.Kernel.Gaussian(kernel.sharpness, kernel.extent, verbose);
% kernelObj = Recon.SysModel.Kernel.Triangle(kernel.extent, verbose);
% kernelObj = Recon.SysModel.Kernel.KaiserBessel(kernel.sharpness, kernel.extent, verbose);
% kernelObj = Recon.SysModel.Kernel.Sinc(kernel.sharpness, kernel.extent, verbose);

% Choose Proximity object
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);
% proxObj = Recon.SysModel.Proximity.L1Proximity(kernelObj, verbose);
clear kernelObj;

% Create System model
systemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
    output_image_size, proxObj, verbose);
% systemObj = Recon.SysModel.ExactSystemModel(traj, overgrid_factor, output_image_size, verbose);

% Choose density compensation function (DCF)
% dcfObj = Recon.DCF.Analytical3dRadial(traj, verbose);
dcfObj = Recon.DCF.Iterative(systemObj, nPipeIter, verbose);
% dcfObj = Recon.DCF.Voronoi(traj, pfile, verbose);
% dcfObj = Recon.DCF.Unity(traj, verbose);

% Choose Reconstruction Model
reconObj = Recon.ReconModel.LSQGridded(systemObj, dcfObj, verbose);
reconObj.crop = 1;
reconObj.deapodize = 1; % Cant deapodize analytical sinc DC = 0

% Reconstruct Data
reconVol = reconObj.reconstruct(pfile.data, traj);

%% Show the result
imslice(abs(reconVol));

%% Create segmentation mask
nClusters = 3;
dims = size(reconVol);
disp('Calculating clusters');
[cluster_idx, cluster_center] = kmeans(abs(reconVol(:)),nClusters,'distance','sqEuclidean', ...
    'Replicates',3);
cluster_idx = reshape(cluster_idx,dims);

% Select lung
h = imslice(cluster_idx);
ax = gca();
disp('Click on each lung (left/right), then press enter...')
[y_idx, x_idx] = getpts(ax);
x_idx = round(x_idx);
y_idx = round(y_idx);
gui = getappdata(h);
slice_dim = gui.UsedByGUIData_m.dimension_selection.Value;
slice_idx = round(gui.UsedByGUIData_m.slice_slider.Value);
im_val = gui.UsedByGUIData_m.hImage.CData;
lung_idx  = im_val(x_idx(1),y_idx(1));

% Find lung by connected components
disp('Finding connected components');
[CC, NUM] = bwlabeln(cluster_idx==lung_idx, 18);
CC_slice = calcImageSlice(Volume(CC), slice_idx, slice_dim);
CC_lung_idx1 = CC_slice(x_idx(1),y_idx(1));
CC_lung_idx2 = CC_slice(x_idx(2),y_idx(2));

% Calculate mask from ventilation volume
[pathstr,name,ext] = fileparts(dixon_pfile);
load([pathstr filesep() name '_gas_recon.mat']);
[ventMask, cluster_center] = kmeans(abs(gasVol(:)),2,'distance','sqEuclidean', ...
    'Replicates',2);
ventMask = reshape(ventMask,dims);
vent_slice = calcImageSlice(Volume(ventMask), slice_idx, slice_dim);
vent_lung_idx1 = vent_slice(x_idx(1),y_idx(1));
vent_lung_idx2 = vent_slice(x_idx(2),y_idx(2));
ventMask = (ventMask(:)==vent_lung_idx1) | (ventMask(:)==vent_lung_idx2);
ventMask = reshape(ventMask,dims);
ventVolume = sum(ventMask(:));

% Create initial mask
lung_mask = (CC == CC_lung_idx1) | (CC == CC_lung_idx2);

% Close lung to not include noise or air outside patient
disp('Performing closing operation...');
lung_closed = zeros([size(lung_mask)]);
NUM = 1;
radius = 0.5;
delta_radius = 0.25;
bigVolumeDifference = true;
while((NUM < 2) | (bigVolumeDifference)) %Perform until the background disappears
    radius = radius + delta_radius;
    
    [xgrid, ygrid, zgrid] = meshgrid(-ceil(radius):ceil(radius));
    ball = (sqrt(xgrid.^2 + ygrid.^2 + zgrid.^2) <= radius);
    
    openedimage = imopen(lung_mask,ball);
    
    
    [CC2, NUM] = bwlabeln(openedimage, 18);
    CC2_slice = calcImageSlice(Volume(CC2), slice_idx, slice_dim);
    CC2_lung_idx1 = CC2_slice(x_idx(1),y_idx(1));
    CC2_lung_idx2 = CC2_slice(x_idx(2),y_idx(2));
    
    % Create mask
    lung_closed = (CC2 == CC2_lung_idx1) | (CC2 == CC2_lung_idx2);
    
    closedVolume = sum(lung_closed(:));
    bigVolumeDifference = abs(1-(ventVolume/closedVolume)) > 0.4;
    
    test = 1;
end
lung_mask2 = lung_mask & lung_closed;

% erode just to make sure we dont get any air phase...
[xgrid, ygrid, zgrid] = meshgrid(-1:1);
ball = (sqrt(xgrid.^2 + ygrid.^2 + zgrid.^2) <= 1);
lung_mask = imerode(lung_mask2,ball);

%% Show the result
imslice(abs(reconVol).*~lung_mask,'Segmented thoracic cavity');

%% save the results
[pathstr,name,ext] = fileparts(bhute_pfile);
uteVol = reconVol;
save([pathstr filesep() name '_bhute_recon.mat'],'uteVol');
niiname = [pathstr filesep() name '_bhute_recon.nii'];
nii = make_nii(abs(reconVol));
save_nii(nii,niiname);

save([pathstr filesep() name '_lungMask.mat'],'lung_mask');
niiname = [pathstr filesep() name '_lungMask.nii'];
nii = make_nii(abs(lung_mask));
save_nii(nii,niiname);
end