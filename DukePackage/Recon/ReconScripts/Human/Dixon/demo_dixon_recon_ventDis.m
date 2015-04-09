
% Required parameters
% 1. Define reconstruction parameters
output_image_size = 128*[1 1 1];
overgrid_factor = 3;
gasKernel.sharpness = 0.33;
gasKernel.extent = 9*gasKernel.sharpness;
dissolvedKernel.sharpness = 0.17;
dissolvedKernel.extent = 9*dissolvedKernel.sharpness;
% kernel.extent = 1.5;
verbose = 1;
nPipeIter = 10;

pfile_path = filepath('C:\Users\Scott\Desktop\subject002_065\P34304.7')

% Human Ventilation Parameters
pfileOverride = GE.Pfile.Pfile();
pfileOverride.rdb.rdb_hdr_user1  = 0.512; % pw_gxwa
pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
pfileOverride.rdb.rdb_hdr_user44 = 1.536; % pw_gxw/1000
pfileOverride.rdb.rdb_hdr_user22 = 0.145; %toff
pfileOverride.rdb.rdb_hdr_user23 = 101; % primeplus
% pfileOverride.rdb.rdb_hdr_user23 = 137.508; % primeplus

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

%% Split pfile into 2: disolved and gas
dissolved_pfile = pfile;
dissolved_pfile.data = dissolved_pfile.data(:,1:2:end);
dissolved_pfile.rdb.rdb_hdr_user20 = size(dissolved_pfile.data,2);

gas_pfile = pfile;
gas_pfile.data = gas_pfile.data(:,2:2:end);
gas_pfile.rdb.rdb_hdr_user20 = size(gas_pfile.data,2);

%% Calculate trajectories (will be same for both pfiles)
% Calculate trajectory for a single radial ray 
radialDistance = MRI.Trajectories.Centric.Radial.calcRadialRay(dissolved_pfile, delays, output_image_size);

% 	% Only keep data during gradients on
[junk, dissolved_pfile] = MRI.DataProcessing.removeNonReadoutSamples(radialDistance, dissolved_pfile);   
[radialDistance, gas_pfile] = MRI.DataProcessing.removeNonReadoutSamples(radialDistance, gas_pfile);   

% Distribute rays onto 3d sphere
traj = MRI.Trajectories.Centric.Distribute.calculate3dTrajectories(radialDistance, dissolved_pfile);

% Undo loopfactor
[junk, dissolved_pfile] = MRI.DataProcessing.undoloopfactor(traj, dissolved_pfile);
[traj, gas_pfile] = MRI.DataProcessing.undoloopfactor(traj, gas_pfile);

% Calculate Maximum volume size for Nyquist
MRI.DataProcessing.calculateNyquistMatrixSize(radialDistance, dissolved_pfile);

% Vectorize data and traj for recon
[junk, dissolved_pfile] = MRI.DataProcessing.vectorizeDataAndTraj(traj, dissolved_pfile);
[traj, gas_pfile] = MRI.DataProcessing.vectorizeDataAndTraj(traj, gas_pfile);

% Enforce Nyquist limits
[junk, dissolved_pfile] = MRI.DataProcessing.enforceNyquistBounds(traj, dissolved_pfile);
[traj, gas_pfile] = MRI.DataProcessing.enforceNyquistBounds(traj, gas_pfile);

%% Reconstruct Gas data
% Choose kernel
gasKernelObj = Recon.SysModel.Kernel.Gaussian(gasKernel.sharpness, gasKernel.extent, verbose);

% Choose Proximity object
gasProxObj = Recon.SysModel.Proximity.L2Proximity(gasKernelObj, verbose);
clear gasKernelObj;

% Create System model
gasSystemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
    output_image_size, gasProxObj, verbose);

% Choose density compensation function (DCF)
gasDcfObj = Recon.DCF.Iterative(gasSystemObj, nPipeIter, verbose);

% Choose Reconstruction Model
gasReconObj = Recon.ReconModel.LSQGridded(gasSystemObj, gasDcfObj, verbose);
gasReconObj.crop = 1;
gasReconObj.deapodize = 1; 

gasVol = gasReconObj.reconstruct(gas_pfile.data, traj);

%% Reconstruct Dissolved Data
dissolvedKernelObj = Recon.SysModel.Kernel.Gaussian(dissolvedKernel.sharpness, dissolvedKernel.extent, verbose);

% Choose Proximity object
dissolvedProxObj = Recon.SysModel.Proximity.L2Proximity(dissolvedKernelObj, verbose);
clear dissolvedKernelObj;

% Create System model
dissolvedystemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
    output_image_size, dissolvedProxObj, verbose);

% Choose density compensation function (DCF)
dissolvedDcfObj = Recon.DCF.Iterative(dissolvedystemObj, nPipeIter, verbose);

% Choose Reconstruction Model
dissolvedReconObj = Recon.ReconModel.LSQGridded(dissolvedystemObj, dissolvedDcfObj, verbose);
dissolvedReconObj.crop = 1;
dissolvedReconObj.deapodize = 1; 

dissolvedVol = dissolvedReconObj.reconstruct(dissolved_pfile.data, traj);

%% Show the results
imslice(abs(gasVol),'Gas');
imslice(abs(dissolvedVol),'Dissolved');

% save the result
[pathstr,name,ext] = fileparts(pfile_path);
save([pathstr filesep() name '_gas_recon.mat'],'gasVol');
save([pathstr filesep() name '_dissolved_recon.mat'],'dissolvedVol');
