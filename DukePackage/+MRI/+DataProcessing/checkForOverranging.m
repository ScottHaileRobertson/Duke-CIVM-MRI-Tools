% CHECKFOROVERRANGING
%
% This function removes remove baseline views from data, and returns an
% accurate pfile.
%
% Note #1: This function assumes that the data is in matrix form [npts x nframes]
%
% Usgae: [pfile] = removeBaselineViews([pfile or pfile_name])
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
function checkForOverranging(varargin)
% Parse inputs
if(nargin < 1)
    [file, path] = uigetfile('*.*', 'Select Pfile');
    pfile_name = strcat(path, file);
    
    % Read pfile header
    pfile = GE.Pfile.Header.read(pfile_name);
else
    if(isa(varargin{1},'GE.Pfile.Pfile'))
        pfile = varargin{1};
    else
        pfile_name = varargin{1};
        
        % Read pfile header
        pfile = GE.Pfile.Header.read(pfile_name);
    end
end

% Check if extended dynamic range is used
switch(pfile.rdb.rdb_hdr_point_size)
	case 2
		% Extended dynamic range is off
		max_value = min(abs(intmax('int16')),abs(intmin('int16')));
	case 4
		% Extended dynamic range is on
		max_value = min(abs(intmax('int32')),abs(intmin('int32')));
	otherwise
		error('Only 2 and 4 are accepted as extended dynamic range options.');
end

if(any(real(pfile.data)>=max_value) | any(imag(pfile.data)>=max_value))
	warning('Data appears to be clipped, you likely overranged!');
end
end