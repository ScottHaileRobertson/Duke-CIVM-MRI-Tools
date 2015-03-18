% This function removes data that was sampled after the decay ramp
%
% Note #1: This function assumes that the data is in matrix form [npts x nframes]
%
% Usage: [correctedData, pfile] = removeNonReadoutSamples(rawData, pfile)
%
% Author: Scott Haile Robertson
% Date: 8/10/2014
%
function [radialDist, pfile] = ...
	removeNonReadoutSamples(radialDist, pfile)
% Pull relavent info out of pfile
bw = pfile.rdb.rdb_hdr_user12;                 % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw);                                             % Time between each sample
dwell_time = Math.nearestMultipleOf(dwell_time,0.002); % Dwell time must be an integer multible of 2us
grad_delay_time = pfile.rdb.rdb_hdr_user22;    % Time between ADC on and start of gradient ramp
ramp_time = pfile.rdb.rdb_hdr_user1;           % Time to ramp gradients
decay_time = pfile.rdb.rdb_hdr_user38;
plat_time = pfile.rdb.rdb_hdr_user44;
npts = pfile.rdb.rdb_hdr_frame_size;

% Put everything in pixel units rather than time
grad_delay_npts = grad_delay_time/dwell_time;
ramp_npts = ramp_time/dwell_time;
plat_npts = plat_time/dwell_time;
decay_npts = decay_time/dwell_time;

% Create sample vector
pts_vec = [0:(npts-1)]';

% Calculate sample number of each region boundary
ramp_start_pt = grad_delay_npts; % Can handle negative gradient start times
plateau_start_pt = ramp_start_pt + ramp_npts;
decay_start_pt = plateau_start_pt + plat_npts;
decay_end_pt = decay_start_pt + decay_npts;

% Check if there are sample points after the decay ramp
start_nonreadout_pts = ceil(decay_end_pt);
if(start_nonreadout_pts < npts)
	% Warn user that we sampled post decay ramp
	warning(['Data was sampled beyone the readout trapezoid... ' ...
		'throwing away ' num2str(npts - start_nonreadout_pts) ' samples. '...
		'Please decrease the number of points to avoid this.']);
	
	% Remove data and radial traj that occur after decay ramp
	radialDist(start_nonreadout_pts:end,:) = [];
	pfile.data(start_nonreadout_pts:end,:)=[];
	
	% Update header to keep it in sync with data
	pfile.rdb.rdb_hdr_frame_size = start_nonreadout_pts-1;
end
end
