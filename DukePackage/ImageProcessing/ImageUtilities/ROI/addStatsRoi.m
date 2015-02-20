%ADDSTATSROI   Adds an ROI that keeps statistics
%
%   addStatsRoi(shape,axes) - Adds a stats roi of the given shape to the
%                             given axes
%
%   addStatsRoi(shape) - Adds a stats roi of the given shape to the current
%                        axes (current axes defined by gca).
%
%   addStatsRoi() - Adds a rectantular stats roi to the current axes
%                   (current axes defined by gca).
%
%   See also: GCA
%
%   Copyright: 2012 Scott Haile Robertson.
%   Website: www.ScottHaileRobertson.com
%   $Revision: 1.0 $  $Date: Dec 10, 2012 $
function roi = addStatsRoi(varargin)
if(nargin > 0)
    % If user specifies an axes, use it.
    ax = varargin{2};
else
    % Default axes is current axes.
    ax = gca();
end

if(nargin > 1)
    % If user specifies a shape, use it.
    switch varargin{1}
        case 'rect'
            roi = imrect(ax);
        case 'ellipse'
            roi = imellipse(ax);
        case 'poly'
            roi = impoly(ax);
        case 'freehand'
            roi = imfreehand(ax);
        otherwise
            error('Shape not supported. Supported shapes: rect, ellipse, poly, freehand.');
    end
else
    % Default shape is a rectangle
    roi = imrect(ax);
end

% Get the image from the axis
im = getimage(ax);

% Add position listener
addNewPositionCallback(roi,@updateROIstats);

% Create Text Object
statsText = text('BackgroundColor',[1 1 1]);

% Initialize the stats
updateROIstats();

    % This function updates the roi stats.
    function [mean_val, stdev_val] = updateROIstats(varargin) 
        mask = roi.createMask;
        roi_data = im(mask);
        roi_data = roi_data(:);
        roi.mean = mean(roi_data);
        roi.std  = std(roi_data);
        roi.min  = min(roi_data);
        roi.max  = max(roi_data);
        roi.area = num2str(sum(mask(:)))

        pos = roi.getPosition();
        set(statsText,'String',['\mu=' sprintf('%4.4f',roi.mean) ' \sigma=' ...
            sprintf('%4.4f',roi.std) ' Range=' sprintf('%4.4f',roi.min) ...
            '-' sprintf('%4.4f',roi.max) ' area=' sprintf('%4.4f',roi.area)]);
        set(statsText,'Position',[pos(1) pos(2)]);
    end
end