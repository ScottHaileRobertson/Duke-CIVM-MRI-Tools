function [ out_vol ] = BF_ND( vol, kernelWidth, r_sig, d_sig )
% N-Dimensional Bilateral Filter (works in parallel)
% Ref:
% C. Tomasi and R. Manduchi
% "Bilateral Filtering for Gray and Color Images"
% Proceedings of the 1998 IEEE International Conference on Computer Vision
%
% Inputs:
% vol: grayscale image to filter
% radius: filter radius in each dimmension
% r_sig: photometric standard deviation in each dimmension
% d_sig: geometric standard deviation in each dimmension
%
% Output:
% vol: filtered version of input "vol"
%
% Author: Scott Haile Robertson (based on 2D version coded by Darin Clark)
% Website: www.ScottHaileRobertson.
%

% Save size before converting to vector
dims = size(vol);
nDim = length(kernelWidth);
sizeDim = 2*kernelWidth+1;
nMesh = prod(sizeDim);

% Initialize matrices
sumD = zeros(size(vol));
out_vol = zeros(size(vol));

% Iterate around filter rather than volume to minimize for loops
dist = zeros([nMesh nDim]);
parfor iDim = 1:nDim
    tmpArray = ones([1 nDim]);
    tmpArray(iDim) = sizeDim(iDim);
    tmpMesh = reshape(-kernelWidth(iDim):kernelWidth(iDim),tmpArray);
    tmpArray = sizeDim;
    tmpArray(iDim) = 1;
    tmpMesh = repmat(tmpMesh,tmpArray);
    dist(:,iDim) = tmpMesh(:);
end
clear tmpArray tmpMesh;

% Precalculate spatial part (does not change based on image)
D_spatial = ones([size(dist,1) 1]);
parfor iDim = 1:nDim
    D_spatial = D_spatial.*exp(-dist(:,iDim)/(2*d_sig(:,iDim)^2));
end

% Precalculate 1/r_sig^2
inv_double_rsig_sqr = 1./(2*r_sig.^2);
for iMesh = 1:nMesh
    disp(['Bilateral filtering' num2str(iMesh) '/' num2str(nMesh)]);

	% Create shifted volume
    shiftVec = zeros([nDim 1]);
    for iDim = 1:nDim
        shiftVec(iDim) = dist(iMesh,iDim);
    end
    shifted_vol = circshift(vol,shiftVec);

	% calculate intensity part
    D_intensity = exp(-single(shifted_vol-vol).^2*inv_double_rsig_sqr); % apply Gaussian kernel
	
    % Accumulate
	out_vol = out_vol + shifted_vol.*D_spatial(iMesh).*D_intensity;
    shifted_vol = [];% Cant use clear in parfor
	sumD = sumD + D_spatial(iMesh)*D_intensity;
end
clear D_spatial D_intensity shifted_vol vol dist;

out_vol = out_vol./sumD;
clear sumD;
end