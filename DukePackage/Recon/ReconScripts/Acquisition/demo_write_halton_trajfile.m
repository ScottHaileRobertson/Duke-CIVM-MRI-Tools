% Works on 2T
% Define the scan parameters
nFrames = 1000000;
filename = ['trajpoints_halton_' num2str(nFrames) '.dat'];

% Use halton sequence to pseudo-randomize
p = haltonset(2);
haltonSeq = net(p,nFrames);

% Calculate spherical coord
azimuthal = 2*pi*haltonSeq(:,1);
elevational = acos(1-2*haltonSeq(:,2));

% Calculate Cartesian coord
x = (cos(azimuthal).*sin(elevational)); %x-coordinates of rays
y = (sin(azimuthal).*sin(elevational)); %y-coordinates of rays
z = (cos(elevational)); %z-coordinates of rays

% Create traj
traj = [x(:)'; y(:)'; z(:)'];
traj = traj(:);

% Write file
fileID = fopen(filename,'w');
fwrite(fileID,traj,'double',0,'ieee-be'); % scanner is BIG ENDIAN
fclose(fileID);
