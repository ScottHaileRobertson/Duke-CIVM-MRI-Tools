function [traj, pfile] = vectorizeDataAndTraj(traj, pfile)
	% Pull relavant info from header
	nDims = 3;
	nPts = pfile.rdb.rdb_hdr_frame_size;
	nFrames = pfile.rdb.rdb_hdr_user20;
	
	% Vectorize data
	pfile.data = reshape(pfile.data,[nPts*nFrames 1]);
	traj = reshape(traj,[nPts*nFrames nDims]);
end