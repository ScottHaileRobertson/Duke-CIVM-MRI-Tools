function demo_dixon_phaseCorrect_splitRes(varargin)

if(nargin < 3 | ~isnumeric(varargin{1}) | ~exist(varargin{2}) | ~exist(varargin{3}))
    RBC_barrier_spect = input('What was the RBC:barrier ratio from spectroscopy?:');
    disp('Select BHUTE pfile');
    bhute_pfile = filepath('/home/scott/Desktop/')
    [pathstr,name,ext] = fileparts(bhute_pfile);
    disp('Select Dixon pfile');
    dixon_pfile = filepath(pathstr);
else
    RBC_barrier_spect = varargin{1};
    bhute_pfile = varargin{2};
    dixon_pfile = varargin{3};
end

[pathstr,name,ext] = fileparts(bhute_pfile);
load([pathstr filesep() name '_lungMask.mat']);
[pathstr,name,ext] = fileparts(dixon_pfile);
load([pathstr filesep() name '_gas_recon.mat']);
load([pathstr filesep() name '_dissolved_recon.mat']);

pfile = GE.Pfile.read(dixon_pfile);
displayPfileHeaderInfo(dixon_pfile);


enforceInLungOnly = 1;

% Create better lung mask from ventilated region.
lungMask = abs(dissolvedVol);
lungMask = lungMask/max(lungMask(:));
lungMask = lungMask>0.1;

radius = 1;
[xgrid, ygrid, zgrid] = meshgrid(-radius:radius);
ball = (sqrt(xgrid.^2 + ygrid.^2 + zgrid.^2) <= radius);
lung_mask = imdilate(lung_mask,ball);
lungMask = lungMask & lung_mask;

%% Finds the phase offset to give the correct spectroscopy ratio
desired_angle = atan2(RBC_barrier_spect,1);
if(enforceInLungOnly)
    netVec = sum(dissolvedVol(lungMask));
else
    netVec = sum(dissolvedVol(:));
end
current_angle_complete = atan2(imag(netVec),real(netVec));
delta_angle_complete = desired_angle-current_angle_complete;

% % Check that angle gives correct RBC:barrier ratio
rotVol = dissolvedVol.*exp(1i*delta_angle_complete);
if(enforceInLungOnly)
    finalVec = sum(rotVol(lungMask));
else
    finalVec = sum(rotVol(:));
end
rbc_barrier_ratio_check1 = imag(finalVec)/real(finalVec)

%% Correct for B0 inhomogeneities
% Itterate until mean phase is zero
gas2 = gasVol;
iterCount = 0;
meanphase = inf;
while((abs(meanphase) > 1E-7))
    if(iterCount > 20)
        warning('Could not get zero mean phase within 10 iterations...');
    end
    if(iterCount > 100)
        error('Could not get zero mean phase within 10 iterations...');
    end
    iterCount = iterCount + 1;
    diffphase = angle(gas2);
    meanphase = mean(diffphase(lungMask(:)));
    gas2 = gas2*exp(-1i*meanphase);
end
diffphase = angle(gas2);

% % remove B0 inhomogeneities
dissolvedVol = dissolvedVol.*exp(1i*(delta_angle_complete-diffphase));

if(enforceInLungOnly)
    finalVec2 = sum(dissolvedVol(lungMask));
else
    finalVec2 = sum(dissolvedVol(:));
end
rbc_barrier_ratio_check2 = imag(finalVec2)/real(finalVec2)

rbc_vol = imag(dissolvedVol);
barrier_vol = real(dissolvedVol);

% Force positiveness and scale 
rbc_vol(rbc_vol<0) = 0;
barrier_vol(barrier_vol<0) = 0;
if(enforceInLungOnly)
    rbc_barrier_ratio_check3 = sum(rbc_vol(lungMask))/sum(barrier_vol(lungMask))
else
    rbc_barrier_ratio_check3 = sum(rbc_vol(:))/sum(barrier_vol(:))
end


maxVal = max(barrier_vol(:));
rbc_vol = rbc_vol/maxVal;
barrier_vol = barrier_vol/maxVal;


% Calculate RBC/Barrier only where there is barrier signal
barrier_thresh = 0.1*max(abs(barrier_vol(:)));
signalMask = abs(barrier_vol)>barrier_thresh;
rbc_to_barrier = zeros(size(barrier_vol));
rbc_to_barrier(lungMask) = rbc_vol(lungMask)./barrier_vol(lungMask);
% rbc_to_barrier(~lungMask)=0;

inLungVox = sum(sum(lungMask,2),3);
rbc_lung = rbc_vol;
rbc_lung(~lungMask) = 0;
barrier_lung = barrier_vol;
barrier_lung(~lungMask) = 0;
rbc_proj = sum(sum(rbc_lung,2),3);
bar_proj = sum(sum(barrier_lung,2),3);
rbc2bar_proj = rbc_proj./bar_proj;
withinLung = inLungVox > 50;

% Should read fov/matrix size from header
dixonPfile = GE.Pfile.read(dixon_pfile);
fov = dixonPfile.rdb.rdb_hdr_fov;
matSize = size(barrier_vol,1);
SI_units = fov*(1:matSize)/matSize;
SI_units_inLung = SI_units(withinLung);
meanSI = mean(SI_units_inLung);
SI_units_inLung = SI_units_inLung - meanSI; % make zero center
SI_units = SI_units - meanSI;
rbc2bar_proj_inLung = rbc2bar_proj(withinLung);

figure();
subplot(2,1,1)
plot(SI_units_inLung,rbc2bar_proj_inLung);
ylabel('RBC:barrier');
xlabel('S/I Position (mm)');
subplot(2,1,2);
plot(SI_units_inLung,inLungVox(withinLung),'.')
ylabel('# Voxels in slice');
xlabel('S/I Position (mm)');

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

projected = struct();
projected.rbc2barrier = rbc2bar_proj_inLung;
projected.si_position = SI_units_inLung;
save([pathstr filesep() name '_projected_SI_Ratios.mat'],'projected');
