%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A demonstration of basic Dynamic Contrast Enhanced (DCE) reconstruction
% using radial data using a simple sliding window reconstruction.
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% code_path='/Users/james/Desktop/DCE_proto';
% run([code_path '/Duke-CIVM-MRI-Tools/setup.m']);
% run([code_path '/GE-MRI-Tools/setup.m']);
% run([code_path '/Non-Cartesian-Reconstruction/setup.m'])
% u_dir='/Users/james/';

%% Slow (but decent) Reconstruction parameters
scale = 1;
oversampling = 1; % Use at least 2. Turn "crop" off on reconObject, then increase this until you dont see signal at the edge of FOV, then turn "crop" back on and use that oversampling
sharpness = 0.3;  % This is a key parameter that tradesoff SNR and resolution (making sharpness smaller will blurr the object, but increase SNR and vice versa)
extent = 9*sharpness; % 9 is a good value

verbose = 0;
nPipeIter = 3; % Use 3-15 for option 1, lots more (~75) for iterative recon. This should scale recon time linearly(ish)
crop = 0;


%% Find directory containing all information
reconDir = ['/home/scott/code/dev/BrukerRecon/B03180/'];

%% Read the header file to get scan info
% This could be done nicer by reading the header, etc. - I was lazy and hard-coded
nPts = 64;
nCoils = 4;
nRaysPerKey = 1980;
nKeys = 13;
nAcq = 11;
samplesPerAcq = nPts*nRaysPerKey*nKeys;
unscaled_size = 2*nPts*[1 1 1];
scaled_output_size = round(scale*unscaled_size);
scale = scaled_output_size(1)/unscaled_size(1);
overgrid_mat_size = scaled_output_size*oversampling;
if(crop)
    reconMatSize = scaled_output_size;
else
    reconMatSize = overgrid_mat_size;
end

%% Sliding window parameters
keysPerWindow = nPts*nRaysPerKey*13; % 10 keys of data (I like 10-15 for this dataset)
windowStep = nPts*nRaysPerKey*1; % Step by one key of data

%% Read in fid data, put in 
dataFile = [ reconDir 'fid'];
fid = fopen(dataFile);
data = fread(fid,inf,'int32');
fclose(fid);
data = complex(data(1:2:end),data(2:2:end)); % Make the data actually complex
data = reshape(data,[nPts nCoils nRaysPerKey nKeys nAcq]); % put data into matrix form
data = permute(data,[1 3 4 5 2]);% Make coil last dimmension
data = reshape(data,[nPts nRaysPerKey*nKeys nAcq nCoils]); % Vectorize all but coil dimmension

%% Read in trajectory
trajFile = [ reconDir  'traj'];
fid = fopen(trajFile);
traj = fread(fid,inf,'double');
fclose(fid);
traj = reshape(traj,[3 nPts nRaysPerKey nKeys]); % put trajectory into matrix form
traj = permute(traj,[2 3 4 1]); % Put [kx,ky,kz] dimmension last
traj = reshape(traj,[nPts nRaysPerKey*nKeys 3])/scale; % vectorize keys and Acq

% Vectorize data and traj
traj = reshape(traj,[nPts*nRaysPerKey*nKeys 3]); % vectorize all but [kx,ky,kz] dimmension
data = reshape(data,[nPts*nRaysPerKey*nKeys*nAcq nCoils]); % Vectorize all but coil dimmension

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
startMod = mod(windowStartIdxs-1,samplesPerAcq)+1;
[uniqueStarts,ia,ic] = unique(startMod); % find all unique trajectories

% and yet he nevr uses either edges or counts.
try % new code 2014b and newer
[uniqueCounts, uniqueEdges] = histcounts(ic,length(uniqueStarts));
catch me % old codeb pre. 
[uniqueCounts, uniqueEdges] = histc(ic,length(uniqueStarts));
end


% Presort traj so parfor works...
store = struct;
nSysMat = length(uniqueStarts);
% store.deltaInput = ones([keysPerWindow 1]);
% deapVol = zeros(overgrid_mat_size);
store.tmpVol = zeros(overgrid_mat_size);
disp(['Completed 0/' num2str(nSysMat) ' traj subsets']);
store.windowTraj = zeros([keysPerWindow 3]);


for iSysMat = 1:nSysMat % This can be done in parallel, but takes LOTS of memory
    % Make a trajectory for each window
    store.windowTraj = squeeze(traj(mod(uniqueStarts(iSysMat)+[0:(keysPerWindow-1)]-1,samplesPerAcq)+1,:));
        
    % Construct system model for this trajectory
    disp('   Creating System model');
    store.systemObj = Recon.SysModel.MatrixSystemModel(store.windowTraj, oversampling, ...
        scaled_output_size, proxObj, verbose);   % This can be stored
    
    % Calculate (Pipe) Iterative Density Compensation weights
    disp('   Calculating DCF');
    store.dcfObj = Recon.DCF.Iterative(store.systemObj, nPipeIter, verbose); % This can be stored
     
%     % Compute deapodization volume for this traj
%     disp('   Calculating Deapodization');
%     deapVol = systemObj'*(store.deltaInput.*dcfObj.dcf);
%     deapVol = reshape(full(deapVol),overgrid_mat_size); % make unsparse;
%     deapVol = ifftshift(ifftn(deapVol));
    
    % Create a data matrix of all repetitions of this trajectory
    sameStartIdx = find(ic==iSysMat); 
    nSameStart = length(sameStartIdx);
    store.dataIdxRep = repmat(windowStartIdxs(sameStartIdx),[keysPerWindow 1]) + repmat([0:(keysPerWindow-1)]',[1 nSameStart]);
    store.windowData = reshape(data(store.dataIdxRep(:),:),[keysPerWindow nSameStart*nCoils]);
    store.dcfRep = repmat(store.dcfObj.dcf,[1 nSameStart*nCoils]);
%     clear dcfObj;
    
    % Grid all data that share this trajectory
    disp(['   Gridding (' num2str(nSameStart) ' time points)x(' num2str(nCoils) ' coil channels) datasets...']);
    store.windowData = store.windowData.*store.dcfRep;
    store.ATrans = store.systemObj.ATrans;
    store.windowRecon = store.ATrans*store.windowData; % We will get a huge speek boost if you can get this to take more advantage of CPU
    
    % Perform SOS across coil channels;
    disp(['   Performing IFFT and SOS on ' num2str(nSameStart) ' time points and ' num2str(nCoils) ' coils...']);
    store.windowRecon = reshape(store.windowRecon, [size(store.systemObj.A,2) nSameStart nCoils]);
%     clear systemObj;

    for iSameStart = 1:nSameStart 
        % Figure out this windows index
        iWindow = sameStartIdx(iSameStart);
        
        for iCoil = 1:nCoils
            % Reconstruct image domain with IFFT
            store.tmpVol = reshape(full(store.windowRecon(:,iSameStart, iCoil)),overgrid_mat_size); % make unsparse;
            store.tmpVol = ifftshift(ifftn(store.tmpVol));
            
            % Accumulate SOS
            slidingWindowReconVol(:,:,:,iWindow) = slidingWindowReconVol(:,:,:,iWindow) + store.tmpVol.^2;%(store.tmpVol.*conj(store.tmpVol));
            disp(['      Finished Coil ' num2str(iCoil) '/' num2str(nCoils)]);
        end
                
        % Finish SOS
        slidingWindowReconVol(:,:,:,iWindow) = sqrt(slidingWindowReconVol(:,:,:,iWindow));
        
        % Show some progress
        disp(['   Completed ' num2str(iSameStart) '/' num2str(nSameStart) ' time Points']); 
    end
     
    % Show some progress
    disp(['Completed ' num2str(iSysMat) '/' num2str(nSysMat) ' traj subsets']);
end
reconTime = toc

% Show the reconstruction
imslice(abs(slidingWindowReconVol),'Sliding window');
