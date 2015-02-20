function writeTiffImage(im, filename) 
% Function to write a 32 bit tiff image from a given image matrix
% 
% Writen by Scott Haile Robertson 05/22/12  
%
% Input    im       = Image matrix to write to file
%          filename = full path to *.out file
%
%%
    %Create Tiff object
    t = Tiff(filename,'w');
    
    %Fill in some info to help with writing the file
    tagstruct.ImageLength = size(im,1);
    tagstruct.ImageWidth = size(im,2);
    tagstruct.Photometric = 1; %Min is Black
    tagstruct.BitsPerSample = 32; %int32
    tagstruct.SampleFormat = Tiff.SampleFormat.Int; % int32
    tagstruct.SamplesPerPixel = 1; %black and white (not RGB)
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    t.setTag(tagstruct);
    
    %write the file
    t.write(int32(im));
    t.close();
end