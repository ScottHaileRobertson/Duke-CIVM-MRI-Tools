% Load volumes
disp('Load Ventilation volume');
% vent_path = filepath('/home/scott/Desktop/')
vent_path = filepath('/home/scott/Public/data/20150514/Subj002_067/P03072_gas_recon.mat')
load(vent_path);

disp('Load Dissolved volume');
% dissolved_path = filepath('/home/scott/Desktop/')
dissolved_path = filepath('/home/scott/Public/data/20150514/Subj002_067/P03072_dissolved_recon.mat')
load(dissolved_path);

%% Correct for B0 inhomogeneities
unwrapped = angle(gasVol);
unwrapped2=unwrap_phase_laplacian(angle(gasVol)); %not using this. 
diffphase = mean(unwrapped2(:))-unwrapped;
diffphase2 = mean(unwrapped(:))-unwrapped;

diss_phasemaped = dissolvedVol.*exp(1i.*diffphase);

% Display
imslice(abs(dissolvedVol),'Pre Bo');
imslice(abs(diss_phasemaped),'Post Bo');

%% Finds the phase offset to give the correct spectroscopy ratio
RBC_barrier_spect = 0.46;

sum_vol = sum(dissolvedVol(:));
sum_vol = sum_vol/abs(sum_vol);
current_angle = atan2(imag(sum_vol),real(sum_vol));
desired_angle = atan2(RBC_barrier_spect,1);
delta_angle = desired_angle-current_angle;

phasedVol = diss_phasemaped*exp(1i*delta_angle);

% Phase sweep

nAng = 360;
angs = linspace(0,2*pi,nAng+1);
angs = angs(1:(end-1));
b = zeros([size(diss_phasemaped) nAng]);
for iAng = 1:nAng
    iAng
    b(:,:,:,iAng) = diss_phasemaped*exp(1i*angs(iAng));
end

imslice(-real(b),'Phase Sweep Real')
imslice(imag(b),'Phase Sweep imag')




