% Get Pfile
pfile_path = filepath('/home/scott/Desktop/');
% pfile_path = filepath('/home/scott/Public/data/20150506/Subj002_066/P11776.7');
% pfile_path = filepath('/home/scott/Desktop/');

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);


%            Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
dis_fit_guess = [   1           -35            215          0; % Component #1
                    1           -285           200          0; % Component #2
                    1           -393           115          0; % Component #3
                    1           -3840           50          0; % Component #4
                    1           -3870           50          0; % Component #4
    ];
% dis_fit_guess = [   1           -37            175          0; % Component #1
%                     1           -360           125          90; % Component #2
%                     1           -3845           50          0; % Component #3
%     ];
% dis_fit_guess = [];

%            Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
gas_fit_guess = [   1           8           15          0]; % Component #

zeropadsize = 2048;
linebroadening = 0; %Hz
skipDownstreamFrames = 25;
minte = 0E-6;
TEs = [875 975 1075 1175]*1E-6;
endToff = 1;
throwAwayCalFrames = 0;
delta_rec_freq = 0;%pfile.image.user12; % freq_off in Hz

% nComp = 3;
% rbc_comp_num = 1;
% barrier_comp_num = 2;
nComp = 5;
rbc_comp_num = 1;
barrier_comp_num = 3;

% Pull relavent info from header
npts = pfile.rdb.rdb_hdr_frame_size;
nFrames = pfile.rdb.rdb_hdr_user20;
bw = pfile.rdb.rdb_hdr_user12;                 % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw*1000);
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
pfile.rdb.rdb_hdr_user12 = 1/(2*dwell_time*1000);
bw = pfile.rdb.rdb_hdr_user12;       
t = (dwell_time*((1:npts) - 1)); %sec

% Remove toff points
t = t(endToff:end);
pfile.data = pfile.data(endToff:end,:);
pfile.rdb.rdb_hdr_frame_size = size(pfile.data,1);
npts = pfile.rdb.rdb_hdr_frame_size;

% Split disolved data into separate pfile
nDisFrames = 200;
nGasFrames = 1;
dis_pfile = pfile;
dis_pfile.data = pfile.data(:,1:nDisFrames);
dis_pfile.rdb.rdb_hdr_user20 = nDisFrames; % nframes
center_freq = pfile.rdb.rdb_hdr_ps_mps_freq/10;

% Split dissolved data into 4 pfiles (one for each TE)
setTE90 = [];%2035E-6;
nTE = 4;
disTE_pfile = cell(1,nTE);
b = zeros(4,2048);
te90 = zeros(nTE,1);
if(nComp > 3)
    plasma_phase = zeros(nTE,1);
end
for iTE = 1:nTE
      
    disTE_pfile{iTE} = dis_pfile;
    disTE_pfile{iTE}.data = disTE_pfile{iTE}.data(1:end,(skipDownstreamFrames*nTE+iTE):nTE:end);
    disTE_pfile{iTE}.rdb.rdb_hdr_frame_size = size(disTE_pfile{iTE}.data,1); %npts
    disTE_pfile{iTE}.rdb.rdb_hdr_user20 = size(disTE_pfile{iTE}.data,2); % nframes
    
    % Undo TE Phase from off center excitation
%     delta_phase = delta_rec_freq*(TEs(iTE)-TEs(1)+t)';
%     delta_phase = repmat(exp(-1i*2*pi*delta_phase),[1 disTE_pfile{iTE}.rdb.rdb_hdr_user20]);
%     disTE_pfile{iTE}.data = disTE_pfile{iTE}.data.*delta_phase;
     delta_phase = delta_rec_freq*(TEs(iTE)-TEs(1));
    delta_phase = exp(-1i*2*pi*delta_phase);
    disTE_pfile{iTE}.data = disTE_pfile{iTE}.data*delta_phase;
        
    % Average
    avgdata = mean(disTE_pfile{iTE}.data,2);
        
    if(isempty(dis_fit_guess))
        nmrMix = NMR_Mix(avgdata, t,[],[],[],[],zeropadsize, linebroadening, center_freq);
        nmrMix = nmrMix.autoAddComponents(nComp);
    else
        nmrMix = NMR_Mix(avgdata, t, dis_fit_guess(:,1),dis_fit_guess(:,2),...
            dis_fit_guess(:,3),dis_fit_guess(:,4),zeropadsize,linebroadening,center_freq);
        nmrMix = nmrMix.fitTimeDomainSignal();
    end
    
    figure();
    nmrMix.plotFit();
    title(['TE' num2str(iTE)]);
    
    nmrMixes{iTE} = nmrMix;
    
    % Calculate TE 90
    deltaF(iTE) = nmrMix.freq(barrier_comp_num)-nmrMix.freq(rbc_comp_num);
    deltaPhase = nmrMix.phase(barrier_comp_num)-nmrMix.phase(rbc_comp_num);
    
    time180(iTE) = abs(1/(2*deltaF(iTE)));
    time90(iTE) = 0.5*time180(iTE);
    
    te90(iTE) = (90-deltaPhase)/(360*deltaF(iTE)) + TEs(iTE);
    
    relTe90 = te90(iTE) - TEs(iTE);
    rbcComp = NMR_Mix([0 0], [0 0],nmrMix.amp(rbc_comp_num),...
        nmrMix.freq(rbc_comp_num),nmrMix.fwhm(rbc_comp_num),...
        nmrMix.phase(rbc_comp_num),zeropadsize, linebroadening, center_freq);
    
    barrierComp = NMR_Mix([0 0], [0 0],nmrMix.amp(barrier_comp_num),...
        nmrMix.freq(barrier_comp_num),nmrMix.fwhm(barrier_comp_num),...
        nmrMix.phase(barrier_comp_num),zeropadsize, linebroadening, center_freq);
    
    if(isempty(setTE90))
        rbc_barrier(iTE) = abs(rbcComp.calcTimeDomainSignal(te90(iTE)))/...
            abs(barrierComp.calcTimeDomainSignal(te90(iTE)));
        totSig = nmrMix.calcTimeDomainSignal(te90(iTE));
        real_imag(iTE) = real(totSig)/imag(totSig);
    else
        rbc_barrier(iTE) = abs(rbcComp.calcTimeDomainSignal(setTE90))/...
            abs(barrierComp.calcTimeDomainSignal(setTE90));
        totSig = nmrMix.calcTimeDomainSignal(setTE90);
        real_imag(iTE) = real(totSig)/imag(totSig);
    end
end

for(iTE = 1:nTE)
    while(te90(iTE)>(minte+time180(iTE)))
        % This te is too high, so subtract 180 deg of phase
        te90(iTE) = te90(iTE) - time180(iTE);
    end
    while(te90(iTE)<minte)
        % This TE is too low, so add 180 deg of phase
        te90(iTE) = te90(iTE) + time180(iTE);
    end
end


% Calculate ratios
colors = get(groot,'defaultAxesColorOrder');
for iTE = 1:nTE
    area_rbc(iTE) = nmrMixes{iTE}.amp(1);
    area_barrier(iTE) = nmrMixes{iTE}.amp(2)+nmrMixes{iTE}.amp(3);
    area_gas(iTE) = nmrMixes{iTE}.amp(4)+nmrMixes{iTE}.amp(5);
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


%% Split Gas spectra into separate pfile
gas_pfile = pfile;
gas_pfile.data = pfile.data(:,nDisFrames + (1:nGasFrames));
gas_pfile.rdb.rdb_hdr_user20 = nGasFrames; % nframes
bw = pfile.rdb.rdb_hdr_user12;
npts = pfile.rdb.rdb_hdr_frame_size;
gasMix = NMR_Mix(gas_pfile.data,t,gas_fit_guess(:,1),gas_fit_guess(:,2),...
    gas_fit_guess(:,3),gas_fit_guess(:,4),zeropadsize,linebroadening,center_freq);
gasMix = gasMix.fitTimeDomainSignal();
% gasMix = NMR_Mix(gas_pfile.data,t,[],[],[],[],zeropadsize,linebroadening,center_freq);
% gasMix = gasMix.fitTool();
figure()
gasMix.plotFit();
gasMix

%% Split flip calibration frames into separate pfile
flipCal_pfile = pfile;
flipCal_pfile.data = pfile.data(:,(nDisFrames+nGasFrames+throwAwayCalFrames+1):end);
nFlipCal = size(flipCal_pfile.data,2);
flipCal_pfile.rdb.rdb_hdr_user20 = nFlipCal; % nframes

% Pick a reasonable DC sample
dc_sample_idx = 5;

% Calculate flip angle
MRI.DataProcessing.calcFlipAngle(flipCal_pfile, dc_sample_idx);

%% Display fits, TE90 etc
te1Mix = nmrMixes{1}
te2Mix = nmrMixes{2}
te3Mix = nmrMixes{3}
te4Mix = nmrMixes{4}

deltaF
mean_te90 = mean(te90);
stdev_te90 = std(te90);
phase90_usec = round(time90*1E6)
te90_usec = round(te90'*1E6)
disp(['TE90=' num2str(mean_te90*1E6) 'usec (' num2str(stdev_te90*1E6) 'usec stdev)']);

if(barrier_comp_num > 2)
    mean_plasma_phase = mean(plasma_phase);
    std_plasma_phase = std(plasma_phase);
    disp(['Plasma component will be ' num2str(mean_plasma_phase) ...
        'degrees (' num2str(std_plasma_phase) ...
        'degree std) out of phase with RBC']);
end

rbc_barrier
disp(['mean RBC:Barrier = ' num2str(mean(abs(rbc_barrier))) ' (' num2str(std(abs(rbc_barrier))) ' std dev)']);
real_imag
disp(['mean Real:Imag = ' num2str(mean(abs(real_imag))) ' (' num2str(std(abs(real_imag))) ' std dev)']);