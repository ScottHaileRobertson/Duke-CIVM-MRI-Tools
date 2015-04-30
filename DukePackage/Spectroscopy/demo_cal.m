pfile_path = filepath('C:\Users\Scott\Desktop\');

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

%            Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
fit_guess = [   1           -20            215          0; % Component #1
                1           -300           200          0; % Component #2
                1           -380           125          0; % Component #3
                1           -3821           50          0; % Component #4
                1           -3845           50          0; % Component #5
];
% fit_guess = [];

zeropadsize = 512;
linebroadening = 0; %Hz
% nComp = 5;
% rbc_comp_num = 1;
% barrier_comp_num = 3;
nComp = 5;
rbc_comp_num = 1;
barrier_comp_num = 3;

% Split disolved data into separate pfile
nDisFrames = 200;
nGasFrames = 1;
dis_pfile = pfile;
dis_pfile.data = pfile.data(:,1:nDisFrames);
dis_pfile.rdb.rdb_hdr_user20 = nDisFrames; % nframes
center_freq = pfile.rdb.rdb_hdr_ps_mps_freq/10;

% Split dissolved data into 4 pfiles (one for each TE)
nTE = 4;
minte = 738;
TEs = [875 975 1075 1175];
disTE_pfile = cell(1,nTE);
b = zeros(4,2048);
te90 = zeros(nTE,1);
if(nComp > 3)
    plasma_phase = zeros(nTE,1);
end
for iTE = 1:nTE
    disTE_pfile{iTE} = dis_pfile;
    disTE_pfile{iTE}.data = disTE_pfile{iTE}.data(:,iTE:nTE:end);
    disTE_pfile{iTE}.rdb.rdb_hdr_user20 = size(disTE_pfile{iTE}.data,2); % nframes
    npts = disTE_pfile{iTE}.rdb.rdb_hdr_frame_size;
    
    % Average
    avgdata = mean(disTE_pfile{iTE}.data,2);
    
    % Pull relavent info from header
    bw = pfile.rdb.rdb_hdr_user12;
    npts = pfile.rdb.rdb_hdr_frame_size;
    
    % 	% Line broaden
    dwell_time = 62E-6;%1/(2*bw*1000);
    broaddata = avgdata;
    t = (dwell_time*((1:npts) - 1)); %sec
    
    if(isempty(fit_guess))
        nmrMix = NMR_Mix(broaddata, t,[],[],[],[],zeropadsize, linebroadening, center_freq);
        nmrMix = nmrMix.autoAddComponents(nComp);
    else
        nmrMix = NMR_Mix(broaddata, t, fit_guess(:,1),fit_guess(:,2),...
            fit_guess(:,3),fit_guess(:,4),zeropadsize,linebroadening,center_freq);
        nmrMix = nmrMix.fitTimeDomainSignal();
    end
    
    figure();
    nmrMix.plotFit();
    title(['TE' num2str(iTE)]);
    
    nmrMixes{iTE} = nmrMix;
    
    % Calculate TE 90
    time180 = 180E6/(360*(nmrMix.freq(rbc_comp_num)-nmrMix.freq(barrier_comp_num)));
    startingPhaseDiff = nmrMix.phase(rbc_comp_num)-nmrMix.phase(barrier_comp_num);
    teDiffFunc = @(t)startingPhaseDiff+360*t*(nmrMix.freq(1)-nmrMix.freq(barrier_comp_num))-90;
    fitoptions = optimoptions('fsolve','Display','off');
    te90(iTE) = fsolve(teDiffFunc,0,fitoptions)*1E6+TEs(iTE);
    deltaPhase = 0;
    while(te90(iTE)<minte)
        te90(iTE) = te90(iTE) + time180;
        deltaPhase = deltaPhase + 180;
    end
    while(te90(iTE)>(minte+time180))
        te90(iTE) = te90(iTE) - time180;
        deltaPhase = deltaPhase - 180;
    end
    if(barrier_comp_num > 2)
        t_te90 = (te90(iTE)-TEs(iTE))*1E-6;
        startingPhaseDiff = nmrMix.phase(rbc_comp_num)-nmrMix.phase(2);
        plasma_phase(iTE) = startingPhaseDiff+360*t_te90*(nmrMix.freq(1)-nmrMix.freq(2))-deltaPhase;
    end
end

te1Mix = nmrMixes{1}
te2Mix = nmrMixes{2}
te3Mix = nmrMixes{3}
te4Mix = nmrMixes{4}

mean_te90 = mean(te90); %usec
stdev_te90 = std(te90);


% Calculate ratios
colors = get(groot,'defaultAxesColorOrder');
% area_rbc = nmrMixes{iTE}.amp(1).*nmrMixes{iTE}.fwhm;
% area_barrier = nmrMixes{iTE}.amp.*nmrMixes{iTE}.fwhm;
% area_gas = nmrMixes{iTE}.amp.*nmrMixes{iTE}.fwhm;
for iTE = 1:nTE
    area_rbc(iTE) = nmrMixes{iTE}.amp(1).*nmrMixes{iTE}.fwhm(1);
    area_barrier(iTE) = nmrMixes{iTE}.amp(2).*nmrMixes{iTE}.fwhm(2);
    area_gas(iTE) = nmrMixes{iTE}.amp(3).*nmrMixes{iTE}.fwhm(3);
end
rbc_barrier = area_rbc./area_barrier;
rbc_gas = area_rbc./area_gas;
gas_barrier = area_gas./area_barrier;
gas_rbc = area_gas./area_rbc;
barrier_rbc = area_barrier./area_rbc;
barrier_gas = area_barrier./area_gas;


figure()
subplot(4,1,1);
plot(repmat(TEs,[3 1])',[area_rbc(:) area_barrier(:) area_gas(:)]);
xlabel('TE');
ylabel('Signal intensity (area under curve)');
legend('rbc','barrier','gas')
subplot(4,1,2);
plot(TEs, rbc_barrier, '-', 'Color',colors(1,:));
hold on;
plot(TEs, gas_barrier, '-', 'Color',colors(3,:));
hold off;
xlabel('TE');
ylabel('Ratio with barrier phase');
legend('rbc:barrier','gas:barrier')
subplot(4,1,3);
plot(TEs, rbc_gas, '-', 'Color',colors(1,:));
hold on;
plot(TEs, barrier_gas, '-', 'Color',colors(2,:));
hold off;
xlabel('TE');
ylabel('Ratio with gas phase');
legend('rbc:gas','barrier:gas');
subplot(4,1,4);
plot(TEs, gas_rbc, '-', 'Color',colors(3,:));
hold on;
plot(TEs, barrier_rbc, '-', 'Color',colors(2,:));
hold off;
xlabel('TE');
ylabel('Ratio with RBC phase');
legend('gas:rbc','barrier:rbc');


% Split Gas spectra into separate pfile
gas_pfile = pfile;
gas_pfile.data = pfile.data(:,nDisFrames + (1:nGasFrames));
gas_pfile.rdb.rdb_hdr_user20 = nGasFrames; % nframes
bw = pfile.rdb.rdb_hdr_user12;
npts = pfile.rdb.rdb_hdr_frame_size;
% %            Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
fit_guess = [   1           8           15          0
    1           10           -55          0]; % Component #
gasMix = NMR_Mix(gas_pfile.data,t,fit_guess(:,1),fit_guess(:,2),...
    fit_guess(:,3),fit_guess(:,4),zeropadsize,linebroadening,center_freq);
gasMix = gasMix.fitTimeDomainSignal();
% gasMix = NMR_Mix(gas_pfile.data,t,[],[],[],[],zeropadsize,linebroadening,center_freq);
% gasMix = gasMix.fitTool();
figure()
gasMix.plotFit();
gasMix

% Split flip calibration frames into separate pfile
throwAwayCalFrames = 5;
flipCal_pfile = pfile;
flipCal_pfile.data = pfile.data(:,(nDisFrames+nGasFrames+throwAwayCalFrames+1):end);
nFlipCal = size(flipCal_pfile.data,2);
flipCal_pfile.rdb.rdb_hdr_user20 = nFlipCal; % nframes

% Pick a reasonable DC sample
dc_sample_idx = 5;

% Calculate flip angle
MRI.DataProcessing.calcFlipAngle(flipCal_pfile, dc_sample_idx);

disp(['TE90=' num2str(mean_te90) 'usec (' num2str(stdev_te90) 'usec stdev)']);

if(barrier_comp_num > 2)
    mean_plasma_phase = mean(plasma_phase);
    std_plasma_phase = std(plasma_phase);
    disp(['Plasma component will be ' num2str(mean_plasma_phase) ...
        'degrees (' num2str(std_plasma_phase) ...
        'degree std) out of phase with RBC']);
end