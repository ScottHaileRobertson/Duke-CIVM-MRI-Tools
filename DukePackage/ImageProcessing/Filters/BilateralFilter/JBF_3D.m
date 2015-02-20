function [ out_vol ] = JBF_3D( vols, radius, r_sig, d_sig )
% 3D Joint Bilateral Filter
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
%

nVols = length(vols);

% Make sure we have the correct number of range weights
% Consider also allowing domain and radii differences per volume
if(length(r_sig) ~= nVols)
	error('You must give and r_sig value for each volume');
end

% Make sure all volumes are the same size
dims = size(vols{1}); 

sumD = cell(size(vols));
out_vol = cell(size(vols));
for iVol=1:nVols
	% Save size before converting to vector
	if(size(vols{iVol})~=dims)
		error('All volumes must have the same size... ');
	end
	
	% Initialize matrices
	sumD{iVol} = zeros(dims);
	out_vol{iVol} = zeros(dims);
end

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
inv_double_rsig_sqr = cell(size(vols));
for iVol=1:nVols
	inv_double_rsig_sqr{iVol} = 1/(2*r_sig{iVol}^2);
end
pctMult = 100/nMesh;
tic
for iMesh = 1:nMesh
	D_intensity = ones(size(vols{1})); %initialize
	for iVol=1:nVols
		% Create shifted volume for each volume
		shifted_vol = circshift(vols{iVol},[x_mesh(iMesh) y_mesh(iMesh) z_mesh(iMesh)]);

		% calculate intensity part
		D_intensity = D_intensity.*exp(-single(shifted_vol-vols{iVol}).^2*inv_double_rsig_sqr{iVol}); % apply Gaussian kernel
	end
	
	for iVol=1:nVols
		% Create shifted volume for each volume
		shifted_vol = circshift(vols{iVol},[x_mesh(iMesh) y_mesh(iMesh) z_mesh(iMesh)]);
		
		out_vol{iVol} = out_vol{iVol} + shifted_vol.*D_spatial(iMesh).*D_intensity;
		sumD{iVol} = sumD{iVol} + D_spatial(iMesh)*D_intensity;
	disp(['Finished ' num2str(iMesh*pctMult) '%, est. ' num2str(toc*(nMesh-iMesh)/(iMesh*60)) ' minutes remaining...']);
	end
end
clear D_spatial D_intensity padded_vol vol x_mesh y_mesh z_mesh;

for iVol=1:nVols
	out_vol{iVol} = out_vol{iVol}./sumD{iVol};
end
