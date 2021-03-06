
% Required parameters
% 1. Define reconstruction parameters
output_image_size = 128*[1 1 1];
overgrid_factor = 3;
kernel.sharpness = 0.33;
kernel.extent = 7*kernel.sharpness;
% kernel.extent = 1.5;
verbose = 1;
nPipeIter = 10;

pfile_path = filepath('/home/scott/Desktop/rohan_20160121')
% pfile_path = filepath('/home/scott/Desktop/')

% Human Ventilation Parameters
pfileOverride = GE.Pfile.Pfile();

pfileOverride.rdb.rdb_hdr_user1  = 0.992; % pw_gxwa
pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
pfileOverride.rdb.rdb_hdr_user44 = 2.976; % pw_gxw/1000
pfileOverride.rdb.rdb_hdr_user22 = 0.05; %toff
% pfileOverride.rdb.rdb_hdr_user23 = 101; % primeplus
pfileOverride.rdb.rdb_hdr_user23 = 137.508; % primeplus
pfileOverride.rdb.rdb_hdr_user32 = 1;

% Gradient delays
delays.x_delay = 0.000;
delays.y_delay = 0.000;
delays.z_delay = 0.000;

%Optional parameters
revision_override = [];  %Optional override if it can't be automatically read from the pfile

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);

% Convert from Pfile format
pfile = convertLegacyPfile(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% Override header values (optional)
displayPfileHeaderInfo(pfile);
pfile = overridePfile(pfile, pfileOverride);
if(verbose)
    % Display key header info
    displayPfileHeaderInfo(pfile);
end

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

% Choose Proximity object
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);
clear kernelObj;

% Create System model
systemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
    output_image_size, proxObj, verbose);

% Choose density compensation function (DCF)
dcfObj = Recon.DCF.Iterative(systemObj, nPipeIter, verbose);

% Choose Reconstruction Model
reconObj = Recon.ReconModel.LSQGridded(systemObj, dcfObj, verbose);
reconObj.crop = 1;
reconObj.deapodize = 1; % Cant deapodize analytical sinc DC = 0
% reconObj = Recon.ReconModel.ConjugateGradient(systemObj, 10, verbose);
% clear modelObj;
% clear dcfObj;

% Reconstruct Data
reconVol = reconObj.reconstruct(pfile.data, traj);

%% Show the result
imslice(abs(reconVol));

% save the result
[pathstr,name,ext] = fileparts(pfile_path);
niiname = [pathstr filesep() name '_recon.nii'];
nii = make_nii(abs(reconVol));
save_nii(nii,niiname,16);