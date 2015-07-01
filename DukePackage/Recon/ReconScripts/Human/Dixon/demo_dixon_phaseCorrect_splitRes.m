
disp('Select BHUTE pfile');
ute_path = filepath('/home/scott/Desktop/');
[pathstr,name,ext] = fileparts(ute_path);
load([pathstr filesep() name '_lungMask.mat']);

disp('Select dixon pfile');
dix_path = filepath('/home/scott/Desktop/');
[pathstr,name,ext] = fileparts(dix_path);
load([pathstr filesep() name '_gas_recon.mat']);
load([pathstr filesep() name '_dissolved_recon.mat']);


%% The calculated RBC:barrier ratio
RBC_barrier_spect =  0.20475;

%% Finds the phase offset to give the correct spectroscopy ratio
desired_angle = atan2(RBC_barrier_spect,1);
netVec = sum(dissolvedVol(lung_mask));
current_angle_complete = atan2(imag(netVec),real(netVec));
delta_angle_complete = desired_angle-current_angle_complete;

% % Check that angle gives correct RBC:barrier ratio
% rotVol = dissolvedVol.*exp(1i*delta_angle_complete);
% finalVec = sum(rotVol(:));
% rbc_barrier_ratio_check = imag(finalVec)/real(finalVec)

%% Correct for B0 inhomogeneities
% Itterate until mean phase is zero
gas2 = gasVol;
iterCount = 0;
meanphase = inf;
while((abs(meanphase) > 1E-7))
    if(iterCount > 10)
        error('Could not get zero mean phase within 10 iterations...');
    end
    iterCount = iterCount + 1
    diffphase = angle(gas2);
    meanphase = mean(diffphase(lung_mask(:)));
    gas2 = gas2*exp(-1i*meanphase);
end
diffphase = angle(gas2);

% % remove B0 inhomogeneities
dissolvedVol = dissolvedVol.*exp(1i*(delta_angle_complete-diffphase));

% Create lung mask from ventilated region.
lungMask = abs(dissolvedVol);
lungMask = lungMask/max(lungMask(:));
lungMask = lungMask>0.2;

% Force positiveness and scale 
rbc_vol = imag(dissolvedVol);
rbc_vol(rbc_vol<0) = 0;
barrier_vol = real(dissolvedVol);
barrier_vol(barrier_vol<0) = 0;
maxVal = max(max(rbc_vol(:)),max(barrier_vol(:)));
rbc_vol = rbc_vol/maxVal;
barrier_vol = barrier_vol/maxVal;

% Calculate RBC/Barrier only where there is barrier signal
barrier_thresh = 0.01*max(barrier_vol(:));
signalMask = barrier_vol>barrier_thresh;
rbc_to_barrier = zeros(size(barrier_vol));
rbc_to_barrier(signalMask) = rbc_vol(signalMask)./barrier_vol(signalMask);
rbc_to_barrier(~lungMask)=0;

imslice(rbc_vol,'RBC')
imslice(barrier_vol,'barrier')
imslice(rbc_to_barrier,'RBC:barrier')

save([pathstr filesep() name '_rbc.mat'],'rbc_vol');
nii = make_nii(rbc_vol);
save_nii(nii,[pathstr filesep() name '_rbc.nii']);

save([pathstr filesep() name '_barrier.mat'],'barrier_vol');
nii = make_nii(barrier_vol);
save_nii(nii,[pathstr filesep() name '_barrier.nii']);

save([pathstr filesep() name '_rbc2barrier.mat'],'rbc_to_barrier');
nii = make_nii(rbc_to_barrier);
save_nii(nii,[pathstr filesep() name '_rbc2barrier.nii']);
