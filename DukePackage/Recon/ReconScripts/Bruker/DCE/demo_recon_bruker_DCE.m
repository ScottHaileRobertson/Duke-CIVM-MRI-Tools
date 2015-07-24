%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A demonstration of basic Dynamic Contrast Enhanced (DCE) reconstruction
% using radial data using a simple sliding window reconstruction.
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reconstruction parameters
output_image_size = 128*[1 1 1];
oversampling = 3; % Use at least 2.5. Turn "crop" off on reconObject, then increase this until you dont see signal at the edge of FOV, then turn "crop" back on and use that oversampling
sharpness = 0.22;  % This is a key parameter that tradesoff SNR and resolution
extent = 9*sharpness; % Never touch this.
verbose = 1;
nPipeIter = 15; % Use 3-15 for option 1, lots more (~75) for iterative recon. This should scale recon time linearly(ish)
crop = 1;
if(crop)
    reconMatSize = output_image_size;
else
    reconMatSize = output_image_size*oversampling;
end

%% Find directory containing all information
% reconDir = uigetdir();
% reconDir = '/home/scott/Desktop/BrukerRecon/bruker_testdat/14';
reconDir = '/home/scott/Desktop/B03180/';

%% Read the header file to get scan info
% This could be done nicer by reading the header, etc. - I was lazy and hard-coded
nPts = 64;
nCoils = 4;
nRaysPerKey = 1980;
nKeys = 13;
nAcq = 11;

%% Sliding window parameters
keysPerWindow = nPts*nRaysPerKey*13; % 10 keys of data (I like 10-15 for this dataset)
windowStep = nPts*nRaysPerKey*5; % Step by one key of data

%% Read in fid data, put in 
% dataFile = filepath('/home/scott/Desktop/BrukerRecon/bruker_testdat/14/fid');
dataFile = filepath('/home/scott/Desktop/BrukerRecon/B03180/fid');
fid = fopen(dataFile);
data = fread(fid,inf,'int32');
fclose(fid);
data = complex(data(1:2:end),data(2:2:end)); % Make the data actually complex
data = reshape(data,[nPts nCoils nRaysPerKey nKeys nAcq]); % put data into matrix form
data = permute(data,[1 3 4 5 2]);% Make coil last dimmension
data = reshape(data,[nPts*nRaysPerKey*nKeys*nAcq nCoils]); % Vectorize all but coil dimmension

%% Read in trajectory
% trajFile = filepath('/home/scott/Desktop/BrukerRecon/bruker_testdat/8/traj');
trajFile = filepath('/home/scott/Desktop/BrukerRecon/B03180/traj');
fid = fopen(trajFile);
traj = fread(fid,inf,'double');
fclose(fid);
traj = reshape(traj,[3 nPts nRaysPerKey nKeys]); % put trajectory into matrix form
traj = permute(traj,[2 3 4 1]); % Put [kx,ky,kz] dimmension last
traj = reshape(traj,[nPts*nRaysPerKey*nKeys 3]); % vectorize all but [kx,ky,kz] dimmension

%% Create the recon objects that dont change 
tic;
% Construct Gridding kernel
kernelObj = Recon.SysModel.Kernel.Gaussian(sharpness, extent, verbose);

% Construct Proximity object
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);

%% Perform sliding window Reconstruction
ptsPerCoil = nPts*nRaysPerKey*nKeys*nAcq; % all keys, acqs
windowStartIdxs = 1:windowStep:(ptsPerCoil-keysPerWindow);
nWindows = length(windowStartIdxs);
slidingWindowReconVol = zeros([reconMatSize nWindows]);

% Batch reconstructions of the same trajectories to save on DCF
% calculations and the computation of system matrices
startMod = mod(windowStartIdxs-1,keysPerWindow)+1;
[uniqueStarts,ia,ic] = unique(startMod); % find all unique trajectories
[uniqueCounts, uniqueEdges] = histcounts(ic,length(uniqueStarts));

nSysMat = length(uniqueStarts);
tmpVol = zeros([reconMatSize max(uniqueCounts(:))]);
disp(['Completed 0/' num2str(nSysMat) ' traj subsets']);
for iSysMat = 1:nSysMat 
    % Make a reconObj for each window
    windowTraj = squeeze(traj(mod(uniqueStarts(iSysMat)+[0:(keysPerWindow-1)]-1,keysPerWindow)+1,:));
    
    % Construct system model
    systemObj = Recon.SysModel.MatrixSystemModel(windowTraj, oversampling, ...
        output_image_size, proxObj, verbose);   
    
    % Option 1: LSQR recon + Pipe DCF (very fast because each coil uses same DCF weights)
    dcfObj = Recon.DCF.Iterative(systemObj, nPipeIter, verbose);
    reconObj = Recon.ReconModel.LSQGridded(systemObj, dcfObj, verbose);
    reconObj.deapodize = 1;
    
    %     % Option 2: Iterative recon (slower as it iterahappen for each
%     % coil
%     reconObj = Recon.ReconModel.ConjugateGradient(systemObj, nPipeIter, verbose);
    
    % Tell recon to crop image to prescribed FOV 
    clear systemObj dcfObj; % free up memory
    reconObj.crop = crop; 
    
    % Reconstruct all data that shares this trajectory
    sameStartIdx = find(ic==iSysMat);
    nSameStart = length(sameStartIdx);
    parfor iSameStart = 1:nSameStart % Has to index 1 upwards... maybe need to use sorted vals...
        iWindow = sameStartIdx(iSameStart);
        sampleIdx = windowStartIdxs(iWindow) + [0:(keysPerWindow-1)];
        
        % Recon each coil in a SOS sense (parallel recon like SENSE/PILS etc
        % would be smarter)
        for iCoil=1:nCoils
            coilData = squeeze(data(sampleIdx,iCoil));
            
            % Compute SOS recon
            tmpVol(:,:,:,iSameStart) = tmpVol(:,:,:,iSameStart) + reconObj.reconstruct(coilData, windowTraj).^2;%;
        end
    end
    clear reconObj windowTraj; % free up memory
    slidingWindowReconVol(:,:,:,sameStartIdx) = tmpVol(:,:,:,1:nSameStart);
    
    % Show some progress
    disp(['Completed ' num2str(iSysMat) '/' num2str(nSysMat) ' traj subsets']);
end

% Take sqrt at end to minimize computation time
clear data traj windowStartIdxs kernelObj proxObj;
slidingWindowReconVol = sqrt(slidingWindowReconVol);
reconTime = toc

% Bilateral filter the data - this can be very slow, but its worth it!
magVol = abs(slidingWindowReconVol); 
clear slidingWindowReconVol;
% Show the reconstruction
imslice(magVol,'Sliding window');

% subVol = magVol(30:100,30:100,60:80,60:70);
% [filt] = BF_ND( subVol, [3 3 3 4], 9E6, [1.25 1.25 1.25 2]);
% close all; imslice(subVol,'Pre filtering'); imslice(abs(filt),'Filtered Vol');

% tic ;
% proxRange = [1.25 1.25 1.25 0.1];
% [filt] = BF_ND( magVol, round(2.2*proxRange), 9E6, proxRange);
% imslice(filt,'Filtered Vol');
% timeToFilt = toc
% 
