


% Required parameters
% 1. Define reconstruction parameters
output_image_size = 64*[1 1 1];
oversampling = 2;
sharpness = 0.3;
extent = 6*sharpness;
% kernel.extent = 1.5;
verbose = 1;
nPipeIter = 5;

pfile_path = filepath('/home/scott/data/demo/human/ventilation/P11776.7_mu_fav_vent')

% Human Ventilation Parameters
pfileOverride = GE.Pfile.Pfile();

pfileOverride.rdb.rdb_hdr_user1  = 0.508; % pw_gxwa
pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
pfileOverride.rdb.rdb_hdr_user44 = 2.048; % pw_gxw/1000
pfileOverride.rdb.rdb_hdr_user22 = 0.13; %toff
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
pfile = GE.Pfile.read(pfile_path);

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

% Pick a reasonable DC sample
dc_sample_idx = MRI.DataProcessing.calcDCsample(max(radialDistance,[],2));

% Throw away all but last DC sample point
pfile.data = pfile.data(dc_sample_idx:end,:);
radialDistance = radialDistance(dc_sample_idx:end,:);
pfile.rdb.rdb_hdr_frame_size = size(pfile.data,1);

% Distribute rays onto 3d sphere
traj = MRI.Trajectories.Centric.Distribute.calculate3dTrajectories(radialDistance, pfile);

% Undo loopfactor
[traj, pfile] = MRI.DataProcessing.undoloopfactor(traj, pfile);

% % Calculate flip angle
% MRI.DataProcessing.calcFlipAngle(pfile, dc_sample_idx);

% Calculate Maximum volume size for Nyquist
MRI.DataProcessing.calculateNyquistMatrixSize(radialDistance, pfile);

% Calculate SNR weights
weights = MRI.DataProcessing.calcSNRWeights(pfile, dc_sample_idx, 1);

% bunch into subsets
nPerSubset = 50;
[nPts nFrames] = size(pfile.data);
subsetStep = 50;
subsetStarts = 1:subsetStep:500;
nSubsets = length(subsetStarts);
fixedData = zeros(nPts*nPerSubset,1,nSubsets);
fixedTraj = zeros(nPts*nPerSubset,3,nSubsets);
for iSubset=1:nSubsets
    subsetFrames = subsetStarts(iSubset)+(0:(nPerSubset-1));
    tmp = pfile.data(:,subsetFrames);
    fixedData(:,1,iSubset) = reshape(tmp,[nPts*nPerSubset 1]);
    
    tmp = reshape(traj(:,subsetFrames,:),[nPts*nPerSubset 3 1]);
    fixedTraj(:,:,iSubset) = tmp;
end

%% Reconstruct data
% Construct Gridding kernel
kernelObj = Recon.SysModel.Kernel.Gaussian(sharpness, extent, verbose);

% Construct Proximity object
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);

% Construct Temporal System Model
systemObj = TemporalSystemModel(fixedTraj, oversampling, output_image_size, proxObj, ...
    verbose);

dcfObj = Recon.DCF.Iterative(systemObj, nPipeIter, verbose);

[nSamples nCoils nTimes] = size(fixedData);
for iTime=1:nTimes
    for iCoil=1:nCoils
        fixedData(:,iCoil,iTime) = fixedData(:,iCoil,iTime).*dcfObj.dcf{iTime};
    end
end
x0 = systemObj' * fixedData;

param.E = systemObj;
param.y = fixedData;
clear pfile traj fixedTraj systemObj fixedData dcfObj weights tmp radialDistance pfileOverride;
param.W = TV_Temp();
param.lambda = 3;
param.display = 1;
param.nite = 8;
x = x0;
for iIter = 1:3
     x = CS_l1(x, param);
end

% Show the reconstruction
imslice(abs(x),'Sliding window');