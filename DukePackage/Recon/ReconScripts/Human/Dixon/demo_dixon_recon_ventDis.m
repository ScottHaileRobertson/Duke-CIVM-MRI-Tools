
%% Required parameters
% Pfile
pfile_path = filepath('/home/scott/Desktop/');

% Prepare to override values
pfileOverride = GE.Pfile.Pfile();

% % For old 64^3
% output_image_sizeg = 64*[1 1 1];
% output_image_sized = 64*[1 1 1];
% overgrid_factor = 3;
% gasKernel.sharpness = 0.35;
% gasKernel.extent = 9*gasKernel.sharpness;
% dissolvedKernel.sharpness = 0.2;
% dissolvedKernel.extent = 9*dissolvedKernel.sharpness;
% verbose = 1;
% nPipeIter = 25;
% pfileOverride.rdb.rdb_hdr_user1  = 0.252; % pw_gxwa
% pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
% pfileOverride.rdb.rdb_hdr_user44 = 1.024; % pw_gxw/1000
% pfileOverride.rdb.rdb_hdr_user22 = 0.125; %toff

% For new 128^3
output_image_sizeg = 64*[1 1 1];
output_image_sized = 64*[1 1 1];
overgrid_factor = 3;
gasKernel.sharpness = 0.35;
gasKernel.extent = 9*gasKernel.sharpness;
dissolvedKernel.sharpness = 0.25;
dissolvedKernel.extent = 9*dissolvedKernel.sharpness;
verbose = 1;
nPipeIter = 25;
pfileOverride.rdb.rdb_hdr_user1  = 0.512; % pw_gxwa
pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
pfileOverride.rdb.rdb_hdr_user44 = 1.536; % pw_gxw/1000
pfileOverride.rdb.rdb_hdr_user22 = 0.125; %toff

deltaf_gas = 0;
deltaf_dissolved = 0; %hz

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
radialDistanceg = MRI.Trajectories.Centric.Radial.calcRadialRay(dissolved_pfile, delays, output_image_sizeg);
radialDistanced = MRI.Trajectories.Centric.Radial.calcRadialRay(dissolved_pfile, delays, output_image_sized);

% % 	% Only keep data during gradients on
% [junk, dissolved_pfile] = MRI.DataProcessing.removeNonReadoutSamples(radialDistanced, dissolved_pfile);   
% [radialDistance, gas_pfile] = MRI.DataProcessing.removeNonReadoutSamples(radialDistanceg, gas_pfile);   

% Distribute rays onto 3d sphere
trajd = MRI.Trajectories.Centric.Distribute.calculate3dTrajectories(radialDistanced, dissolved_pfile);
trajg = MRI.Trajectories.Centric.Distribute.calculate3dTrajectories(radialDistanceg, gas_pfile);

% Undo loopfactor
[trajd, dissolved_pfile] = MRI.DataProcessing.undoloopfactor(trajd, dissolved_pfile);
[trajg, gas_pfile] = MRI.DataProcessing.undoloopfactor(trajg, gas_pfile);

% Demodulate signal
bw = pfile.rdb.rdb_hdr_user12;                 % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw*1000);                                             % Time between each sample
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
nPts = pfile.rdb.rdb_hdr_frame_size;           % Number of sample points per frame/ray
tMat = repmat(dwell_time*(1:nPts)',[1 size(dissolved_pfile.data,2)]);
dissolved_pfile.data = dissolved_pfile.data.*exp(1i*2*pi*deltaf_dissolved*tMat);
gas_pfile.data = gas_pfile.data.*exp(1i*2*pi*deltaf_gas*tMat);

% Calculate Maximum volume size for Nyquist
MRI.DataProcessing.calculateNyquistMatrixSize(radialDistanced, dissolved_pfile);
MRI.DataProcessing.calculateNyquistMatrixSize(radialDistanceg, gas_pfile);

% Vectorize data and traj for recon
[trajd, dissolved_pfile] = MRI.DataProcessing.vectorizeDataAndTraj(trajd, dissolved_pfile);
[trajg, gas_pfile] = MRI.DataProcessing.vectorizeDataAndTraj(trajg, gas_pfile);

% Enforce Nyquist limits
[trajd, dissolved_pfile] = MRI.DataProcessing.enforceNyquistBounds(trajd, dissolved_pfile);
[trajg, gas_pfile] = MRI.DataProcessing.enforceNyquistBounds(trajg, gas_pfile);

%% Reconstruct Gas data
% Choose kernel
gasKernelObj = Recon.SysModel.Kernel.Gaussian(gasKernel.sharpness, gasKernel.extent, verbose);

% Choose Proximity object
gasProxObj = Recon.SysModel.Proximity.L2Proximity(gasKernelObj, verbose);
clear gasKernelObj;

% Create System model
gasSystemObj = Recon.SysModel.MatrixSystemModel(trajg, overgrid_factor, ...
    output_image_sizeg, gasProxObj, verbose);

% Choose density compensation function (DCF)
gasDcfObj = Recon.DCF.Iterative(gasSystemObj, nPipeIter, verbose);

% Choose Reconstruction Model
gasReconObj = Recon.ReconModel.LSQGridded(gasSystemObj, gasDcfObj, verbose);
gasReconObj.crop = 1;
gasReconObj.deapodize = 1; 

gasVol = gasReconObj.reconstruct(gas_pfile.data, trajg);

%% Reconstruct Dissolved Data
dissolvedKernelObj = Recon.SysModel.Kernel.Gaussian(dissolvedKernel.sharpness, dissolvedKernel.extent, verbose);

% Choose Proximity object
dissolvedProxObj = Recon.SysModel.Proximity.L2Proximity(dissolvedKernelObj, verbose);
clear dissolvedKernelObj;

% Create System model
dissolvedystemObj = Recon.SysModel.MatrixSystemModel(trajd, overgrid_factor, ...
    output_image_sized, dissolvedProxObj, verbose);

% Choose density compensation function (DCF)
dissolvedDcfObj = Recon.DCF.Iterative(dissolvedystemObj, nPipeIter, verbose);1/128

% Choose Reconstruction Model
dissolvedReconObj = Recon.ReconModel.LSQGridded(dissolvedystemObj, dissolvedDcfObj, verbose);
dissolvedReconObj.crop = 1;
dissolvedReconObj.deapodize = 1; 

dissolvedVol = dissolvedReconObj.reconstruct(dissolved_pfile.data, trajd);

%% Show the results
imslice(abs(gasVol),'Gas');
imslice(abs(dissolvedVol),['Dissolved - ' num2str(deltaf_dissolved)]);
% imshowpair(squeeze(abs(gasVol(:,:,52))),squeeze(abs(dissolvedVol(:,:,52))))
title(num2str(deltaf_dissolved))

% save the result
[pathstr,name,ext] = fileparts(pfile_path);
save([pathstr filesep() name '_gas_recon.mat'],'gasVol');
save([pathstr filesep() name '_dissolved_recon.mat'],'dissolvedVol');

[pathstr,name,ext] = fileparts(pfile_path);
niiname = [pathstr filesep() name '_gas_recon.nii'];
nii = make_nii(abs(gasVol));
save_nii(nii,niiname);

[pathstr,name,ext] = fileparts(pfile_path);
niiname = [pathstr filesep() name '_dissolved_recon.nii'];
nii = make_nii(abs(dissolvedVol));
save_nii(nii,niiname);
