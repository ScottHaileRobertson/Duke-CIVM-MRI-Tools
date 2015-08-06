%CALC_ARCHIMEDIAN_SPIRAL_TRAJECTORIES   Archimedian spiral trajectory generator.
%   CALC_ARCHIMEDIAN_SPIRAL_TRAJECTORIES(nframes, primeplus, rad) calculates
%   the radial trajectories following an archimedial spiral. nframes defines
%   the number of frames, primeplus defines the randomizing allong the 
%   spiral, and rad defines the radial sample locations. This is a 
%   slightly modified version of a vectorized version provided by 
%   Sam Johnston.
%
%   Authors: Gary Cofer, Sam Johnston, Scott Haile Robertson.
%   $Revision: 1.0 $  $Date: 2012/07/19 $
function traj = calcArchimedianSpiralTrajectories(radialDistance, header)
% Read relavent info from header
nframes = header.rdb.rdb_hdr_user20;
primeplus = header.rdb.rdb_hdr_user23;
pi_float = single(3.14159265358979323846);

% nframes=floor(nframes/2)*2+1; %Must have odd number of frames
cview=floor(nframes/2)+1;     %Center frame 

is = 0:nframes-1; %In Gary's code, i=acview_start

z_coords = single(abs(1-(is/cview))); %In Gary's code f=fThing
angs = single(primeplus.*is.*(pi_float/180)); %azimuthal Angle in radians
ds = single(sqrt(1-(z_coords.^2)));
x_coords = ds.*single(cos(angs));
y_coords = ds.*single(sin(angs));

%Handle negatives
z_coords = single(z_coords.*((2*(is<=cview))-1));

%normalize
ivec_lengths = single(1./sqrt((x_coords.^2) + (y_coords.^2) + (z_coords.^2)));
xs = single(x_coords.*ivec_lengths);
ys = single(y_coords.*ivec_lengths);
zs = single(z_coords.*ivec_lengths);

thetas = single(acos(zs));
phis = single(atan2(y_coords,x_coords));

dx = single(single(sin(thetas)).*single(cos(phis)));
dy = single(single(sin(thetas)).*single(sin(phis)));
dz = single(single(cos(thetas)));

traj = zeros([size(radialDistance,1) nframes size(radialDistance,2)]);
traj(:,:,1) = radialDistance(:,1)*dx; %x-coordinates of rays
traj(:,:,2) = radialDistance(:,2)*dy; %y-coordinates of rays
traj(:,:,3) = radialDistance(:,3)*dz; %z-coordinates of rays
end