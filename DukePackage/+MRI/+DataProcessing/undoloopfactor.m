function [traj, pfile] = undoloopfactor(traj, pfile)
% Get relavent info from header
loop_factor = pfile.rdb.rdb_hdr_user10;
nframes = pfile.rdb.rdb_hdr_user20;
per_nufft = pfile.rdb.rdb_hdr_user32;

if((per_nufft ~= 1) & (loop_factor > 1))
	% Warn user that loopfactor and non archimedial traj dont work
	h = warning(['If per_nufft~=1, loopfactor cant be greater than 1! ' ...
		'Setting loopfactor to 1...'],'!! Warning !!');
    pfile.rdb.rdb_hdr_user10 = 1;
end

old_idx = 1:nframes;
new_idx = mod((old_idx-1)*loop_factor,nframes)+1;

pfile.data(:,old_idx) = pfile.data(:,new_idx);
traj(:,old_idx, :) = traj(:,new_idx,:);

% Update header in case you try to undo loopfactor again (it will do
% nothing)
pfile.rdb.rdb_hdr_user10 = 1;