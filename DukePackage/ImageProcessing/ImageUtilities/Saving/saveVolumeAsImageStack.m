function saveVolumeAsImageStack(vol, basename, format)
numSlices = size(vol,3);

% Calculate volume min and max
min_val = min(vol(:));
max_val = max(vol(:));

for i=1:numSlices
    file_name = [basename '_' num2str(i) '.' format];
    
    im_matrix = squeeze(vol(:,:,i));
    
    saveMatrixAsImage(im_matrix, file_name, min_val, max_val);
end