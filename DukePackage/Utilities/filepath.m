%FILEPATH   Easy way to get the full path to a file.
%
%   Copyright: 2012 Scott Haile Robertson.
%   Website: www.ScottHaileRobertson.com
%   $Revision: 1.0 $  $Date: Dec 12, 2012 $
function file_path = filepath(file_path, varargin)
% Be nice and ask for file if none is given
if((nargin < 1) || isempty(file_path))
    if(nargin > 1)
        [file, path] = uigetfile('*.*', varargin{1});
    else
        [file, path] = uigetfile('*.*', 'Select file');
    end
    if(isnumeric(file) & isnumeric(path) & file==0 & path ==0)
        file_path = [];
    else
        file_path = strcat(path, file);
    end
end

% Check that the file or directory exists
if(exist(file_path))
    % If its a directory, start looking from there and ask for a file
    if(isdir(file_path))
        [file, path] = uigetfile('*.*', 'Select file',file_path);
        if(isnumeric(file) & isnumeric(path) & file==0 & path ==0)
            file_path = [];
        else
            file_path = strcat(path, file);
        end
    end
else
    file_path = [];
end