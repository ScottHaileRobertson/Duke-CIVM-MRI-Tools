% This function calculates the matrix size that just meets the Nyquist
% limit of the given radial trajectory. 
% Usage: header = calculateNyquistMatrixSize(radialDistance, [pfile])
%
% Author: Scott Haile Robertson
% Date: 8/10/2014
%
function calculateNyquistMatrixSize(radialDistance, pfile)
	% Use the original npts, not after any data processing because this is
	% how the gradient amplitudes are scaled
	npts_orig = pfile.rdb.rdb_hdr_user7;
	nyquistMatrixSize = 4*npts_orig*max(radialDistance(:));
	
	disp(['Matrix Size for Nyquist limit = ' num2str(floor(nyquistMatrixSize))]);
end