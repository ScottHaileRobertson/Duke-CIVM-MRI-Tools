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

% nframes=floor(nframes/2)*2+1; %Must have odd number of frames
cview=floor(nframes/2)+1;     %Center frame 

is = 0:nframes-1; %In Gary's code, i=acview_start

z_coords = abs(1-(is/cview)); %In Gary's code f=fThing
angs = primeplus.*is.*(pi/180); %azimuthal Angle in radians
ds = sqrt(1-(z_coords.^2));
x_coords = ds.*cos(angs);
y_coords = ds.*sin(angs);

%Handle negatives
z_coords = z_coords.*((2*(is<=cview))-1);

%normalize
ivec_lengths = 1./sqrt((x_coords.^2) + (y_coords.^2) + (z_coords.^2));
xs = x_coords.*ivec_lengths;
ys = y_coords.*ivec_lengths;
zs = z_coords.*ivec_lengths;

thetas = acos(zs);
phis = atan2(y_coords,x_coords);

dx = (sin(thetas).*cos(phis));
dy = (sin(thetas).*sin(phis));
dz = (cos(thetas));

traj = zeros([size(radialDistance,1) nframes size(radialDistance,2)]);
traj(:,:,1) = radialDistance(:,1)*dx; %x-coordinates of rays
traj(:,:,2) = radialDistance(:,2)*dy; %y-coordinates of rays
traj(:,:,3) = radialDistance(:,3)*dz; %z-coordinates of rays
end