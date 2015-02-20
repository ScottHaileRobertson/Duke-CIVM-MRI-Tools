% Input parameters
movieName = 'testMovie';
im_size = 100;
nIter = 15;

% Prepare the new file.
vidObj = VideoWriter(movieName);
vidObj.FrameRate = 5;
open(vidObj);

% Set up initial image
fig = figure();
imageHandle = imagesc(zeros(im_size,im_size));
colormap(gray);
axis image;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
    
% Loop and create frames
for i=1:nIter  
    % Update image
    set(imageHandle ,'CData',rand(im_size,im_size));
    
    % Write each frame to the file.
    writeVideo(vidObj,getframe(fig));
end

% Close the file.
close(vidObj);
