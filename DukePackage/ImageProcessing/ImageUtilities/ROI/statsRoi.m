%STATSROI   An ROI that keeps and optionally displays statistics
%
%   statsRoi(shape,axes,...) - Adds a stats roi of the given shape to the
%                             given axes
%
%   statsRoi(shape,...) - Adds a stats roi of the given shape to the current
%                        axes (current axes defined by gca).
%
%   statsRoi(...) - Adds a rectantular stats roi to the current axes
%                   (current axes defined by gca).
%
%   See also: GCA
%
%   Copyright: 2012 Scott Haile Robertson.
%   Website: www.ScottHaileRobertson.com
%   $Revision: 1.0 $  $Date: Dec 13, 2012 $
classdef statsRoi < imroi
    
    properties (SetAccess = 'protected',GetAccess = 'public')
        mean
        std
        min
        max
        area
    end
    
    properties (SetAccess = 'protected',GetAccess = 'protected')
        roi
        text
    end
    
    methods
        function obj = statsRoi(varargin)
            % statsRoi Constructor for statsRoi
            
            % Initialize using super class
            if(nargin > 0)
                %get type of roi
                switch varargin{1}
                    case 'rect'
                        roi = imrect(varargin{2:end});
                    case 'ellipse'
                        roi = imellipse(varargin{2:end});
                    case 'poly'
                        roi = impoly(varargin{2:end});
                    case 'freehand'
                        roi = imfreehand(varargin{2:end});
                    otherwise
                        error('Shape not supported. Supported shapes: rect, ellipse, poly, freehand.');
                end
            else
                roi = imrect(gca());
            end
            
            % Create container object
            obj = obj@imroi(roi.h_group,roi.draw_api);
            obj.roi = roi;
            
            % Add position listener
            addNewPositionCallback(obj,@obj.updateROIstats);
            
            % Create Text Object
            obj.text = text('BackgroundColor',[1 1 1], ...
                'VerticalAlignment','bottom', ...
                'HorizontalAlignment','Left');
            
            % Initialize the stats
            obj.updateROIstats(obj);
            
            % Add to context menu
            contextMenu = obj.roi.api.getContextMenu();
            uimenu(contextMenu,'Label','Show Stats','Callback',@obj.toggleStats);
            uimenu(contextMenu,'Label','Delete','Callback',@obj.handleDeleteEvent);
            obj.roi.api.setContextMenu(contextMenu);
        end
        
        % This function updates the roi stats.
        function [mean_val, stdev_val] = updateROIstats(varargin)
            obj = varargin{1};
            mask = obj.roi.createMask;
            im = getimage(ancestor(obj.roi.h_group,'axes'));
            roi_data = im(mask);
            roi_data = roi_data(:);
            obj.mean = mean(roi_data);
            obj.std  = std(roi_data);
            obj.min  = min(roi_data);
            obj.max  = max(roi_data);
            obj.area = sum(mask(:));
            
            % Update the stats text
            obj.updateStatsText();
        end
        
        function updateStatsText(varargin)
            obj = varargin{1};
            pos = obj.roi.getPosition();
            
            set(obj.text,'String',['\mu=' sprintf('%4.4f',obj.mean) ' \sigma=' ...
                sprintf('%4.4f',obj.std) ' min=' sprintf('%4.4f',obj.min) ...
                ' max=' sprintf('%4.4f',obj.max) ' area=' sprintf('%4.4f',obj.area)]);
            set(obj.text,'Position',[pos(1,1) pos(1,2)]);
        end
        
        function toggleStats(obj, varargin)
            switch get(obj.text,'Visible')
                case 'on'
                    set(obj.text,'Visible','off');
                case 'off'
                    set(obj.text,'Visible','on');
            end
        end
        
        function handleDeleteEvent(obj, varargin)
            obj.roi.delete(); % Let super class handle it
        end
        
        function delete(obj)
            obj.roi.delete(); % Let super class handle it
            delete(obj.text);
        end
        
        function setPosition(obj,pos)
            obj.roi.setPosition(pos); % Let super class handle it
        end
        
        function pos = getPosition(obj)
            pos = obj.roi.getPosition(); % Let super class handle it
        end
        
        function setResizable(obj,TF)
            obj.roi.setResizable(TF); % Let super class handle it
        end
        
        function setFixedAspectRatioMode(obj,TF)
            obj.setFixedAspectRatioMode(TF); % Let super class handle it
        end
        
    end
end

