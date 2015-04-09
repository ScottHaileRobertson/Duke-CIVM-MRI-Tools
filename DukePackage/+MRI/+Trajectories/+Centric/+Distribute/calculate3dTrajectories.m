% This function will appropriately distribute nFrames of radial rays on the
% surface of a sphere. It returns a matrix of trajectories that is 
% [nPts nFrames nDim] and is scaled from 0 to 0.5.
%
% Note #1: This function assumes that the data is in matrix form [npts x nframes]
% Note #2: This function assumes that the traj is in matrix form [npts x nframes x nDims]
%
% Usgae: traj = calculate3dTrajectories(radialDistance, header)
%
% Author: Scott Haile Robertson
% Date: 8/10/2014
%
function traj = calculate3dTrajectories(radialDistance, pfile)
    % Pull relavent info from header
	per_nufft = pfile.rdb.rdb_hdr_user32;
	
    % Choose trajectory calculations based on per_nufft
    if ( per_nufft == 1)
         traj = MRI.Trajectories.Centric.Distribute.ArchimedianSpiral(radialDistance, pfile);
    end
    if (per_nufft == 0)
        traj = double(MRI.Trajectories.Centric.Distribute.GoldenMeans(radialDistance, pfile));
	end
end