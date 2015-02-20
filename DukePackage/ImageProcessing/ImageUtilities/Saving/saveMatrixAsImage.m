function saveMatrixAsImage(im_matrix, file_name, min_val, max_val)
im_matrix = squeeze(im_matrix);
dims = size(im_matrix);

ndims = length(dims);
if(ndims==2)
    im_matrix = repmat(im_matrix,[1,1,3]);
    dims = size(im_matrix);
    ndims = length(dims);
end

if(ndims~=3)
    error('This function only works on 2-D matrices or 3D matrices with RGB channels'); 
end

if(isempty(min_val) & isempty(max_val))
    max_val = max(im_matrix(:));
    min_val = min(im_matrix(:));
end

% Scale image to be between 0 and 1
im_matrix = (im_matrix-min_val)/max_val;

%Save png image file
imwrite(im_matrix,file_name);