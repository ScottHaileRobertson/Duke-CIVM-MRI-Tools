
% Ask user for all pfiles
disp('Locate Phase Calibration Pfile please...');
phaseCal_pfile = filepath();
[pathstr,name,ext] = fileparts(phaseCal_pfile);
disp('Locate Dixon Pfile please...');
dixon_pfile = filepath(pathstr);
disp('Locate BHUTE Pfile please...');
bhute_pfile = filepath(pathstr);

% Perform phase calibration
meanRbc2barrier = demo_dixon_phase_calibration(phaseCal_pfile);


% Reconstruct ventilation/dissolved phase image
demo_dixon_recon_ventDis(dixon_pfile);

% Reconstruct BHUTE and make lung mask 
demo_dixon_recon_BHUTE(bhute_pfile, dixon_pfile)

% Perform dixon decomposition
demo_dixon_phaseCorrect_splitRes(meanRbc2barrier, bhute_pfile, dixon_pfile)
