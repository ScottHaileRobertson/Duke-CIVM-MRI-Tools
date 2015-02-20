%LOADIMAGEDATASET   Loads a set of images into a volume.
%
% Works for dicom images as well as other formats (jpeg, png, etc)
%
%   [vol] = loadImageDataset();
%
%   [vol] = loadImageDataset(image_file_list);
%
%   Copyright: 2012 Scott Haile Robertson.
%   Website: www.ScottHaileRobertson.com
%   $Revision: 1.0 $  $Date: Dec 7, 2012 $
function vol = loadImageDataset(varargin)

% Get image files
if(nargin < 1)
    % If no files are provided, use uigetfile to find them
    [image_files,p_name,filt_idx] = uigetfile( '*.*', 'Select Images', 'MultiSelect', 'on');
    
    % Make image_files the full path name
    numFiles = length(image_files); 
    for i=1:numFiles
        image_files{i} = [p_name image_files{i}];
    end
else
    image_files = varargin{1};
end

% Get file format
[pathstr, name, ext] = fileparts(image_files{1});
isDicom = strcmp(ext,'.dcm');

% Read in volume
num_imgs = size(image_files,2);
if(isDicom)
    % Lookup file info
    im_info = dicominfo(image_files{1});
    
    % Preallocate volume
    vol = zeros(im_info.Width, im_info.Height, num_imgs);
    
    % Load volume
    for i=1:num_imgs
        vol(:,:,i) = dicomread(image_files{i});
    end
else
    % Lookup file info
    im_info = imfinfo(image_files{1});
    
    % Preallocate volume
    vol = zeros(im_info.Width, im_info.Height, num_imgs);
    
    % Load volume
    for i=1:num_imgs
        vol(:,:,i) = imread(image_files{i});
    end
end
end %function