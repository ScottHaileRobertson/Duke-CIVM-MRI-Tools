% SHOWORTHOGONALPLANES  - Displays orthogonal 2D planes of a 3D volume.
%
% Shows triplanar surface renderings, which can each be moved dynamically
% by dragging a plane.
%
% Copyright: 2012 Scott Haile Robertson.
% Website: www.ScottHaileRobertson.com
%
% $Revision: 1.0 $
% $Date: Nov 28, 2012 $
function hImage = showTriplanar(vol, varargin)

% Get initial slices if they are given
if(nargin == 4)
end

% Create initial slices
[x,y,z] = meshgrid(1:size(vol,2),1:size(vol,1),1:size(vol,3));
% slices=slice(x,y,z,vol,size(vol,2)/2,size(vol,1)/2,size(vol,3)/2);
slices=slice(x,y,z,vol,59,60,50);
slice_x=slices(1);
slice_y=slices(2);
slice_z=slices(3);
set(slice_x,'EdgeColor','None','Tag','SliceX');
set(slice_y,'EdgeColor','None','Tag','SliceY');
set(slice_z,'EdgeColor','None','Tag','SliceZ');
colormap(gray);

% Initialize dragging data
clickedAxes = 0;
clickedAxesDim = 0;
active_slice = 0;
start_pt = 0;
cur_data = 0;

% Label axes
xlabel('x');
ylabel('y');
zlabel('z');
xlim([1 size(vol,1)]);
ylim([1 size(vol,2)]);
zlim([1 size(vol,3)]);

% Listen for drag events
addSliceListener(slice_x);
addSliceListener(slice_y);
addSliceListener(slice_z);

    % Sets up the ability for when a slice is being dragged
    function addSliceListener(hSlice)
        set(hSlice,'ButtonDownFcn',@slicePressed);
    end

    % Reacts to the initial click on a slice and prepares for dragging
    function slicePressed(src,evnt)
        thisfig = gcbf();
        set(thisfig,'WindowButtonMotionFcn',@sliceDragged);
        set(thisfig,'WindowButtonUpFcn',@sliceReleased);
        set(thisfig,'ButtonDownFcn','');

        % Get the slice data and dimmension of the clicked slice
        active_slice = src;
        switch get(active_slice,'Tag')
            case 'SliceX'
                cur_data = get(src,'XData');
                clickedAxesDim = 1 %MATLAB calls x the 2nd dim
            case 'SliceY'
                cur_data = get(src,'YData');
                clickedAxesDim = 2;%MATLAB calls x the 1st dim
            case 'SliceZ'
                cur_data = get(src,'ZData');
                clickedAxesDim = 3;
            otherwise
                error('Tag not supported');
        end
        
        % Get the 3D position of the starting click point for drag calcs
        clickedAxes = gca;
        start_pt = get(clickedAxes,'CurrentPoint');
        
        % Update axis being dragged
        updateAxisLabels(round(max(cur_data(:))));
    end

    % Drags a slice as the user moves the mouse
    function sliceDragged(src,evnt)
        % Calculate how much the mouse has moved since the mouse was
        % pressed. There is a slight bug here with the points relative to
        % two different planes.
        cur_pt = get(clickedAxes,'CurrentPoint');
        delta = calcPositionChange(start_pt,cur_pt,clickedAxesDim);
        
        % Calculate the new slice data and the slice index
        newSlice = cur_data+delta;
        new_idx = round(max(newSlice(:)));
        if((new_idx>size(vol,clickedAxesDim)) || (new_idx<1))
            % We're out of bounds, so just stop
            return;
        end
        
        % Update axes and image with new slice
        updateAxisLabels(new_idx);
        newIm = calcImageSlice(vol, new_idx, clickedAxesDim);
        switch clickedAxesDim
            case 1
                set(active_slice,'XData',newSlice);
                set(active_slice,'CData',newIm);
            case 2
                set(active_slice,'YData',newSlice);
                set(active_slice,'CData',newIm);
            case 3
                set(active_slice,'ZData',newSlice);
                set(active_slice,'CData',newIm);
            otherwise
                error('Incalid selected axis dimension.');
        end
    end

    % Stop keeping track of slice when the mouse is released and reset the
    % axes to just x,y,z
    function sliceReleased(src,evnt)
        thisfig = gcbf();
        set(thisfig,'WindowButtonUpFcn','');
        set(thisfig,'WindowButtonMotionFcn','');
        xlabel('x=');
        ylabel('y');
        zlabel('z');
    end

    % Handles some of the weirdness of the 'CurrentPoint' in the axis
    function pt = calcPositionChange(pt1,pt2,dim)
        max_pos = 30;
        if((pt2(2,dim) <= 1) || (pt2(2,dim) >= max_pos))
            pt = pt2(1,dim)-pt1(1,dim);
        else
            pt = pt2(2,dim)-pt1(2,dim);
        end
    end

    % Updates the axis of the currently dragged slice to show the displayed
    % slice number
    function updateAxisLabels(slice)
        switch clickedAxesDim
            case 1
                xlabel(['x=' num2str(slice)]);
                ylabel('y');
                zlabel('z');
            case 2
                xlabel('x');
                ylabel(['y=' num2str(slice)]);
                zlabel('z');
            case 3
                xlabel('x');
                ylabel('y');
                zlabel(['z=' num2str(slice)]);
            otherwise
                error('Dimension not supported');
        end
    end
end
