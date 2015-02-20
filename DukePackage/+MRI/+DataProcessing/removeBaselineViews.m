% This function removes remove baseline views from data, and returns an
% accurate header.
%
% Note #1: This function assumes that the data is in matrix form [npts x nframes]
%
% Usgae: [baselinelessData, header] = removeBaselineViews(rawData, header)
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
function [pfile] = removeBaselineViews(pfile)
	% Pull relavent info out of header
	baseline_skip_size = pfile.rdb.rdb_hdr_da_yres; 
	nframes = pfile.rdb.rdb_hdr_user20;
	npts = pfile.rdb.rdb_hdr_frame_size;
	
	% Remove baselines views
	pfile.data(:, 1:baseline_skip_size:nframes) = []; 

	% Update pfile header to keep it in sync with data
	pfile.rdb.rdb_hdr_user20 = length(pfile.data(:))/npts;
	pfile.rdb.rdb_hdr_da_yres = inf; % In case you try to remove baselines again, nothing will happen
end
