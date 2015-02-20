% Make image of random data
imagesc(rand(100,100));

% add a stats roi. You will need to click on the image to create the roi.
statsRoi('ellipse');

%Try dragging around the roi - the stats update!