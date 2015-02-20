function [ out_vol ] = BF_3D( vol, radius, r_sig, d_sig )
% 3D Bilateral Filter (works in parallel)
% Ref:
% C. Tomasi and R. Manduchi
% "Bilateral Filtering for Gray and Color Images"
% Proceedings of the 1998 IEEE International Conference on Computer Vision
%
% Inputs:
% vol: grayscale image to filter
% radius: filter radius
% r_sig: photometric standard deviation
% d_sig: geometric standard deviation
%
% Output:
% vol: filtered version of input "vol"
%
% Author: Scott Haile Robertson
% Based on 2D version coded by Darin Clark
% Works in parallel 
%

% Save size before converting to vector
dims = int16(size(vol)); % save on memory

% Initialize matrices
sumD = zeros(size(vol));
out_vol = zeros(size(vol));

% Iterate around filter rather than volume to minimize for loops
lsp = int16(-radius:radius);
[x_mesh y_mesh z_mesh] = meshgrid(lsp,lsp,lsp);
clear lsp
r_mesh_sqr = double(x_mesh.^2 + y_mesh.^2 + z_mesh.^2);
within_radius = (r_mesh_sqr <= (radius^2));
x_mesh = x_mesh(within_radius);
y_mesh = y_mesh(within_radius);
z_mesh = z_mesh(within_radius);
clear within_radius 
nMesh = length(x_mesh(:));

% Precalculate spatial part (does not change based on image)
D_spatial = single(exp(-r_mesh_sqr/(2*d_sig^2)));
clear r_mesh_sqr;

% Precalculate 1/r_sig^2
inv_double_rsig_sqr = 1/(2*r_sig^2);
pctMult = 100/nMesh;
tic
parfor iMesh = 1:nMesh
	% Create shifted volume
	shifted_vol = circshift(vol,[x_mesh(iMesh) y_mesh(iMesh) z_mesh(iMesh)]);

	% calculate intensity part
	D_intensity = exp(-single(shifted_vol-vol).^2*inv_double_rsig_sqr); % apply Gaussian kernel
	
	out_vol = out_vol + shifted_vol.*D_spatial(iMesh).*D_intensity;
% 	clear shifted_vol;
	sumD = sumD + D_spatial(iMesh)*D_intensity;
% 	clear D_intensity;
% 	disp(['Finished ' num2str(iMesh*pctMult) '%, est. ' num2str(toc*(nMesh-iMesh)/(iMesh*60)) ' minutes remaining...']);
end
clear D_spatial D_intensity padded_vol vol x_mesh y_mesh z_mesh;

out_vol = out_vol./sumD;
