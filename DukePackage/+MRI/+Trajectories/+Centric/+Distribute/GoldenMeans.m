function traj = calcGoldenMeanTrajectories(radialDistance, header)
% Get relavant info from header 
nFrames = header.rdb.rdb_hdr_user20;

% Has to be single to match up with EPIC float value...
pi_float = single(pi);

i = single(0:nFrames-1);

% Can be calculated with golden_means = single(calc_golden_means(2));
phi1 = single(0.465571224689483642578125);
phi2 = single(0.682327806949615478515625);

tempa= single(mod(i.*phi2,1));
tempb= single(mod(i.*phi1,1));

alpha = single(2*pi_float.* tempa);
beta = single(acos(2*tempb-1)); % Between pos/neg 1

dx = single(cos(alpha).*sin(beta));
dy = single(sin(alpha).*sin(beta));
dz = single(cos(beta));

traj(:,:,1) = radialDistance(:,1)*dx;
traj(:,:,2) = radialDistance(:,2)*dy;
traj(:,:,3) = radialDistance(:,3)*dz;
end