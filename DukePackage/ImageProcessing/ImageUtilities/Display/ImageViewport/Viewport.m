%Viewport This object is used to display one or more volumes on a simgle
%port.
%
%   Copyright: 2012 Scott Haile Robertson.
%   Website: www.ScottHaileRobertson.com
%   $Revision: 1.0 $  $Date: 2013/04/20 $
classdef Viewport < handle
    properties (SetAccess = protected)
        Parent;
        Volumes = []; % Arrray of Volumes
        Colormaps = []; % Array of colormaps for each volume
        Position = []; % Current Volume slice position
        Linkedports = [];
        Rotation = [];
        
        % SubPanels
        TopPanel;
        ImagePanel;
        Image;
        ImageAxis;
        SliderPanel;
        CursorStatsText;
        Sliders = [];
        
        % Layout Defaults
        TopPanelHeight = 14;
        TextSize = 12;
        TopPanelColor = [25, 35, 95]/255;
        TextColor = 'white';
        
        % Listener variables
        MousePressedCursorPosition;
        MousePressedAxisPosition;
        MousePressedParentPosition;
        MousePressedCLim;
    end
    
    methods
        function this = Viewport(varargin)
            % First variable is the Volume(s)
            if(nargin > 0)
                if(isa(varargin{1},'Volume'))
                    this.Volumes = [varargin{1}];
                    varargstart = 2;
                elseif(isempty(varargin{1}))
                    % You dont have to provide a volume
                else
                    this.Volumes = [Volume(varargin{1})];
                    varargstart = 2;
                end
            end
            
            for i=varargstart:2:nargin
                if(isstr(varargin{i}))
                    switch(lower(varargin{i}))
                        case 'parent'
                            this.Parent = varargin{i+1};
                        otherwise
                            error(['Argument ' varargin{i} ' not found']);
                    end
                else
                    error('Arguments are all key-value pairs...');
                end
            end
            
            if(isempty(this.Parent ))
                % If no parent is given, make a figure and use that.
                this.Parent = figure();
            end
            
            % Set up listeners to parent
            set(this.Parent, 'ResizeFcn', @this.resize, ...
                'WindowButtonDownFcn', @this.buttonPressed,...
                'WindowScrollWheelFcn' , @this.mouseScrolled,...
                'WindowButtonUpFcn', @this.buttonReleased,...
                'WindowButtonMotionFcn', @this.pixelHovered,...
                'KeyPressFcn'          , @this.keyPressed);
            %                 'WindowButtonDownFcn'  , @this.buttonPressed, ...
            %                 'WindowButtonMotionFcn', @this.hover, ...
            %
            %             'CloseRequestFcn'      , @closeGUI, ...
            
            
            
            % Set up subpanels
            this.TopPanel = uipanel('BorderType','none',...
                'Units','pixels',...
                'BackgroundColor',this.TopPanelColor,...
                'Parent',this.Parent);
            this.ImagePanel = uipanel('BorderType','none',...
                'Units','pixels',...
                'BackgroundColor','black',...
                'Parent',this.Parent);
            this.ImageAxis = axes('Parent',this.ImagePanel,...
                'XTickLabel',[],...
                'YTickLabel',[],...
                'XTick',[],...
                'YTick',[],...
                'xcolor','green',...
                'ycolor','green',...
                'CLimMode','manual');
            this.CursorStatsText = uicontrol(this.TopPanel, 'Style', 'text', ...
                'String', '', 'Units','pixels', 'FontUnits', 'pixels', ...
                'FontSize',this.TextSize, 'HorizontalAlignment', 'right',...
                'BackgroundColor', this.TopPanelColor,'ForegroundColor',this.TextColor);
                       
            %             initialize position and rotation order
            this.Rotation = 1:this.Volumes(1).NDims;
            this.Position = ones(1,this.Volumes(1).NDims);
            
            % Load Volume if one exists
            if(length(this.Volumes) > 0)
                % Check volume dimensions, add scroll panel if necessary
                if(this.Volumes(1).NDims > 2)
                    this.SliderPanel = uipanel('BorderType','none',...
                        'Units','pixels',...
                        'BackgroundColor','red',...
                        'Parent',this.Parent);
                    
                    for i=3:(this.Volumes(1).NDims)
                        sliderNum = i-2;
                        this.Sliders(sliderNum) = uicontrol(this.SliderPanel, 'Style', 'slider', ...
                            'String', ['Dim' num2str(i)], 'Callback', @slider_callback, ...
                            'Units', 'pixels');
                        addlistener(this.Sliders(sliderNum),'Value','PostSet',@slider_callback);
                        
                        max_v = this.Volumes(1).Dims(i);
                        step_v = [1 min(round(max_v/4),20)]/(max_v-1);
                        
                        set(this.Sliders(sliderNum),'Max',max_v,'Min',1,'Value',1,'SliderStep',step_v);
                    end
                end
            end
            
            
            % Force image update
            this.updateImage();
            this.resize([],[]);
            colormap(gray);
            
            % Initialize WW/WL to scale with volume
            min_v = min(this.Volumes(1).Data(:));
            max_v = max(this.Volumes(1).Data(:));
            if(min_v == max_v)
                if(max_v > 0)
                    min_v = 0;
                else
                    max_v = 1;
                end
            end
            set(this.ImageAxis,'CLim',[min_v max_v],...
                'XTickLabel',[],...
                'YTickLabel',[],...
                'XTick',[],...
                'YTick',[],...
                'xcolor','green',...
                'ycolor','green');
%             cb = double(scribe.colorbar(this.ImageAxis,'EastOutside', [],...
%                 'xcolor','green','ycolor','green'));
%             set(cb,'Visible','off');
%             get(cb)
            
            function slider_callback(hObject,eventdata)
                newPosition = this.Position;
                for i=1:length(this.Sliders)
                    dim = find(this.Rotation == i+2);
                    newPosition(dim) = get(this.Sliders(i),'Value');
                end
                
                this.setPosition(round(newPosition));
            end
        end
        
        function resize(this, eventdata, varargin)
            parent_sz = get(this.Parent, 'Position');
            
            %Build panel layout from bottom left to top right
            %[distFromLeft DistFromBottom Width Height]
            bottomPanelHeight = 0;
            sliderHeight = 20;
            totalSliderPanelHeight = length(this.Sliders)*sliderHeight;
            
            % Add sliders at bottom
            set(this.SliderPanel, 'Position', [1,1,parent_sz(3),length(this.Sliders)*sliderHeight]);
            
            % Add slider bars for extra dims
            for i=1:length(this.Sliders)
                %                     set(h.dim_val_extra_txt(i), 'Position', [1,scroll_panel_height-(i+1)*single_scroll_height,position_width,single_scroll_height]);
                set(this.Sliders(i), 'Position', [1,sliderHeight*(i-1)+1,parent_sz(3),sliderHeight]);
                %                     set(h.dim_name_txt(i), 'Position', [2+scroll_sz(3)-dimmension_sel_width,scroll_panel_height-(i+1)*single_scroll_height,dimmension_sel_width,single_scroll_height]);
            end
            
            
            set(this.ImagePanel, 'Position',[1,...
                bottomPanelHeight+1+totalSliderPanelHeight,...
                parent_sz(3),...
                parent_sz(4)-bottomPanelHeight-this.TopPanelHeight-totalSliderPanelHeight]);
            
            set(this.TopPanel, 'Position',[1,...
                parent_sz(4)-this.TopPanelHeight,...
                parent_sz(3),...
                this.TopPanelHeight]);
            
            % Add stats text
            topPanelSz = get(this.TopPanel, 'Position');
            set(this.CursorStatsText, 'Position',...
                [0 0 topPanelSz(3) (2+this.TextSize)]);
            
        end
        
        function keyPressed(this, eventdata, varargin)
            if(nargin > 0)
                switch varargin{1}.Key
                    case 'uparrow'
                        % Should use rotation matrix
                        oldx = this.Rotation == 1;
                        newx = this.Rotation == 3;
                        this.Rotation(oldx)=3;
                        this.Rotation(newx)=1;
                        this.updateImage();
                    case 'downarrow'
                        % Should use rotation matrix
                        oldx = this.Rotation == 1;
                        newx = this.Rotation == 3;
                        this.Rotation(oldx)=3;
                        this.Rotation(newx)=1;
                        this.updateImage();
                    case 'leftarrow'
                        % Should use rotation matrix
                        oldx = this.Rotation == 2;
                        newx = this.Rotation == 3;
                        this.Rotation(oldx)=3;
                        this.Rotation(newx)=2;
                        this.updateImage();
                    case 'rightarrow'
                        % Should use rotation matrix
                        oldx = this.Rotation == 2;
                        newx = this.Rotation == 3;
                        this.Rotation(oldx)=3;
                        this.Rotation(newx)=2;
                        this.updateImage();
                    otherwise
                        theKey = varargin{1}.Key
                end
            end
        end
        
        function closeGUI(hObject, eventdata, varargin)
        end
        
        function pixelHovered(this, hObject, eventdata, varargin)
            updateCursorStats(this);
        end
        
        function updateCursorStats(this)
            % Get hovered position
            pix = round(get(this.ImageAxis, 'CurrentPoint'));
            x = pix(1,1);
            y = pix(1,2);
            clear pix;
            xdim = this.Volumes(1).Dims(find(this.Rotation == 2));
            ydim = this.Volumes(1).Dims(find(this.Rotation == 1));
            
            %Make sure we're hovering over image
            if(x > 0 & y>0 & x<= xdim & y<=ydim)
                % Get pixel value
                eval(['pixelVal = ' createVolumeIdxString(this, 1, num2str(y), num2str(x)) ';']);
                
                % Make pixel string
                pix_string = ['(' num2str(y) ',' num2str(x)];
                if(this.Volumes(1).NDims>2)
                    for i=3:this.Volumes(1).NDims
                        pix_string = [pix_string ',' num2str(this.Position(find(this.Rotation == i)))];
                    end
                end
                pix_string = [pix_string ')=' num2str(pixelVal) ' '];
                
                % Stop showing stats text
                set(this.CursorStatsText, 'String', pix_string);
            else
                % Stop showing stats text
                set(this.CursorStatsText, 'String', '');
            end
        end
        
        function buttonPressed(this, hObject, eventdata, varargin)
            % Get Start Position of mouse
            this.MousePressedCursorPosition = get(this.Parent,'CurrentPoint');
            this.MousePressedAxisPosition = get(this.ImageAxis,'Position');
            this.MousePressedParentPosition = get(this.Parent,'Position');
            this.MousePressedCLim = get(this.ImageAxis,'CLim');
            
            % Listen for dragging
            switch get(this.Parent, 'SelectionType')
                case 'normal'
                    set(this.Parent, 'WindowButtonMotionFcn', @this.imageDragged);
                case 'alt'
                    set(this.Parent, 'WindowButtonMotionFcn', @this.zoomed);
                case 'extend'
                    set(this.Parent, 'WindowButtonMotionFcn', @this.contrastChanged);
            end
        end
        
        function buttonReleased(this, hObject, eventdata, varargin)
            set(this.Parent, 'WindowButtonMotionFcn' ,@this.pixelHovered);
        end
        
        function zoomed(this, hObject, eventdata, varargin)
            zoom_scale = 500;
            diff_pos = get(this.Parent, 'CurrentPoint') - this.MousePressedCursorPosition;
            zoom_factor = 1+2*(diff_pos(1,1)/zoom_scale);
            if(zoom_factor <= 0)
                return;
            end
            
            new_corners = 0.5-(zoom_factor*(0.5-this.MousePressedAxisPosition));
            set(this.ImageAxis, 'Position', [new_corners(1), ...
                new_corners(2), ...
                zoom_factor*this.MousePressedAxisPosition(3), ...
                zoom_factor*this.MousePressedAxisPosition(4)]);
        end
        
        function contrastChanged(this, hObject, eventdata, varargin)
            contrast_scale = 200;
            diff_pct = (get(this.Parent, 'CurrentPoint') - this.MousePressedCursorPosition)/contrast_scale;
            
            ww = (this.MousePressedCLim(2)-this.MousePressedCLim(1));
            center_val = 0.5*(this.MousePressedCLim(2)+this.MousePressedCLim(1));
            
            % Window Level
            new_center = center_val + diff_pct(1);
            
            % Window Width - Note, cant be zero or negative!
            new_ww = max(ww + diff_pct(2),eps);
            
            % calculate CLims
            new_min = new_center - 0.5*new_ww;
            new_max = new_center + 0.5*new_ww;
            
            set(this.ImageAxis,'CLim',[new_min new_max]);
        end
        
        function imageDragged(this, hObject, eventdata, varargin)
            diff_pos = get(this.Parent, 'CurrentPoint') - this.MousePressedCursorPosition;
            
            newPos = [this.MousePressedAxisPosition(1)+diff_pos(1,1)/this.MousePressedParentPosition(3) ...
                this.MousePressedAxisPosition(2)+diff_pos(1,2)/this.MousePressedParentPosition(4)...
                this.MousePressedAxisPosition(3) ...
                this.MousePressedAxisPosition(4)];
            
            set(this.ImageAxis, 'Position', newPos);
        end
        
        function hover(hObject, eventdata, varargin)
        end
        
        function mouseScrolled(this, hObject, eventdata, varargin)
            % If its a 3+ dimensional volume, scroll through the third
            % dimension
            if(this.Volumes(1).NDims>2)
                this.setPosition(this.Position + eventdata.VerticalScrollCount*(this.Rotation==3));
            end
        end
        
        function setPosition(this, newPosition);
            % Make sure new position has a valid number of dimensions
            if(length(newPosition) ~= this.Volumes(1).NDims)
                error('Must set a position within the volume dimensions');
            end
            
            % Make sure new position is not old position
            if(~isempty(this.Position) & (newPosition == this.Position))
                return;
            end
            
            % Check that each dimmension is within the maximum size of
            % the Volume
            if(any((newPosition < 1) | (newPosition > this.Volumes(1).Dims)))
                return;
            end
            
            % Set the new position
            this.Position = newPosition;
            
            % Update Sliders if necessary
            for i=3:this.Volumes(1).NDims
                newSliderVal = this.Position(this.Rotation == i);
                if(get(this.Sliders(i-2),'Value') ~= newSliderVal )
                    set(this.Sliders(i-2),'Value',newSliderVal );
                end
            end
            
            % Update any listeners
            
            % Update the image
            this.updateImage();
        end
        
        function idxString = createVolumeIdxString(this, volNum, dim1str, dim2str)
            % Create index string
            dim1 = find(this.Rotation == 1);
            dim2 = find(this.Rotation == 2);
            
            % Build string to index slice
            idxString = ['this.Volumes(' num2str(volNum) ').Data('];
            for i=1:length(this.Rotation)
                if(i~=1)
                    idxString = [idxString ','];
                end
                
                if(i == dim1)
                    idxString = [idxString dim1str];
                elseif(i == dim2)
                    idxString = [idxString dim2str];
                else
                    idxString = [idxString num2str(this.Position(i))];
                end
            end
            idxString = [idxString ')'];
        end
        
        function updateImage(this, dim1str, dim2str)
            % Build string to index slice
            idxString = createVolumeIdxString(this, 1, ':', ':');
            idxString = ['slice_im = ' idxString ';'];
            eval(idxString);
            
            if(isempty(this.Image))
                this.Image = imagesc(squeeze(permute(slice_im,this.Rotation)));
            else
                set(this.Image,'CData',squeeze(permute(slice_im,this.Rotation)));
            end
            
            this.updateCursorStats();
            
            axis image;
        end
    end
end