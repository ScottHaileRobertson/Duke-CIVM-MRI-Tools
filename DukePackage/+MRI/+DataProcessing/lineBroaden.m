% This function adds exponential line broadening to the data
%
% Note #1: This function assumes that the data is in matrix form [npts x nframes]
%
% Usage: [pfileWithLineBroadening] = lineBroaden(pfileWithoutLineBroadening, 
%                                                lineBroadeningInHz)
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
function [pfile] = lineBroaden(pfile,lineBroadening)
	% Pull relavent info out of header
	baseline_skip_size = pfile.rdb.rdb_hdr_da_yres; 
	nframes = pfile.rdb.rdb_hdr_user20;
	npts = pfile.rdb.rdb_hdr_frame_size;
	bw = 1000*pfile.rdb.rdb_hdr_user12;                    % Receiver bandwidth (Hz)
    
	% Calculate sample times
    dwell_time = 1/(2*bw);                                 % Time between each sample
    dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
    t = repmat(dwell_time*(0:(npts-1))',[1 nframes]);

	% Create linebroadening vector
    lineBroadeningVec = exp(-lineBroadening*t);
    
	% Line broaden data
    pfile.data = pfile.data.*lineBroadeningVec;
end
