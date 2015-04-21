
%   Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
dis_fit_guess = [
    1           -35            215          0; % Component #1
    1           -285           200          0; % Component #2
    1           -393           115          0; % Component #3
    1           -3840           50          0; % Component #4
    1           -3870           50          0; % Component #4
    ];

% amp_lb = zeros(1,size(dis_fit_guess,1));
amp_lb = [0.05 0.05 0.2 0 0];
freq_lb = [-80 -340 -430 -3890 -4290]
fwhm_lb = [130 130 100 0 0];
phase_lb = -inf*ones(1,size(dis_fit_guess,1));

amp_ub = inf*ones(1,size(dis_fit_guess,1));
freq_ub = [10 -220 -350 -3000 -3820];
fwhm_ub = [270 270 180 35 35];
phase_ub = inf*ones(1,size(dis_fit_guess,1));

gas_fit_guess = [   1           8           15          0]; % Component #

% Miscelaneous parameters
rbc_idx = 1;
barrier_idx = 2:3;
rbc_te90_idx = 1;
barrier_te90_idx = 3;
gas_idx = [4:5];
endToff = 1;
zeropadsize = 2048;
linebroadening = 0; %Hz
skipDownstreamFrames = 25;
throwAwayCalFrames = 0;
minte = 300E-6;
TEs = [875 975 1075 1175]*1E-6;
nDisFrames = 200;
nGasFrames = 1;
nTE = length(TEs);
nComp = size(dis_fit_guess,1);
delta_rec_freq = 0;% pfile.image.user12; % freq_off in Hz

% Get Pfile
pfile_path = filepath('/home/scott/Public/data/')

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

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
t = t(endToff:end)-t(endToff);
pfile.data = pfile.data(endToff:end,:);
pfile.rdb.rdb_hdr_frame_size = size(pfile.data,1);
npts = pfile.rdb.rdb_hdr_frame_size;

%% Split Gas spectra into separate pfile
gas_pfile = pfile;
gas_pfile.data = pfile.data(:,nDisFrames + (1:nGasFrames));
gas_pfile.rdb.rdb_hdr_user20 = nGasFrames; % nframes
bw = pfile.rdb.rdb_hdr_user12;
npts = pfile.rdb.rdb_hdr_frame_size;
gasFit = NMR_Fit(gas_pfile.data,t,zeropadsize,linebroadening,gas_fit_guess(:,1),gas_fit_guess(:,2),...
    gas_fit_guess(:,3),gas_fit_guess(:,4));
gasFit = gasFit.fitTimeDomainSignal();
% gasMix = NMR_Mix(gas_pfile.data,t,[],[],[],[],zeropadsize,linebroadening,center_freq);
% gasMix = gasMix.fitTool();
figure()
gasFit.plotFit();

%% Split flip calibration frames into separate pfile
flipCal_pfile = pfile;
flipCal_pfile.data = pfile.data(:,(nDisFrames+nGasFrames+throwAwayCalFrames+1):end);
nFlipCal = size(flipCal_pfile.data,2);
flipCal_pfile.rdb.rdb_hdr_user20 = nFlipCal; % nframes

% Pick a reasonable DC sample
dc_sample_idx = 5;

% Calculate flip angle
[flip_angle, flip_err] = MRI.DataProcessing.calcFlipAngle(flipCal_pfile, dc_sample_idx);

%% Split disolved data into separate pfile
dis_pfile = pfile;
dis_pfile.data = pfile.data(:,1:nDisFrames);
dis_pfile.rdb.rdb_hdr_user20 = nDisFrames; % nframes
center_freq = pfile.rdb.rdb_hdr_ps_mps_freq/10;

% Consider each TE separately
amplitudes = zeros(nTE,nComp);
frequencies = zeros(nTE,nComp);
fwhms = zeros(nTE,nComp);
phases = zeros(nTE,nComp);
barrier_ratio = zeros(nTE,nComp);
gas_ratio = zeros(nTE,nComp);
ded_gas_ratio = zeros(nTE,nComp);

b = zeros(4,2048);
te90 = zeros(nTE,1);
if(nComp > 3)
    plasma_phase = zeros(nTE,1);
end
for iTE = 1:nTE
    %
    teData = dis_pfile.data(:,(skipDownstreamFrames*nTE+iTE):nTE:end);
    
    % Undo TE Phase from off center excitation
    delta_phase = delta_rec_freq*(TEs(iTE)-TEs(1));
    delta_phase = exp(-1i*2*pi*delta_phase);
    teData = mean(teData,2);
    teData = teData.*delta_phase;
    
    % Calculate fit
    nmrFit = NMR_Fit(teData, t, zeropadsize,linebroadening, dis_fit_guess(:,1),dis_fit_guess(:,2),...
        dis_fit_guess(:,3),dis_fit_guess(:,4));
%     nmrFit = nmrFit.setBounds( amp_lb, amp_ub, freq_lb, freq_ub,...
%         fwhm_lb, fwhm_ub, phase_lb, phase_ub);
    nmrFit = nmrFit.fitTimeDomainSignal();
    
    % save fits
    amplitudes(iTE,:) = nmrFit.nmrMix.amp;
    frequencies(iTE,:) = nmrFit.nmrMix.freq;
    fwhms(iTE,:) = nmrFit.nmrMix.fwhm;
    phases(iTE,:) = nmrFit.nmrMix.phase;
    
    % Calculate ratios
    barrier_ratio(iTE,:)=amplitudes(iTE,:)/sum(amplitudes(iTE,barrier_idx));
    gas_ratio(iTE,:)=amplitudes(iTE,:)/sum(amplitudes(iTE,gas_idx));
    ded_gas_ratio(iTE,:)=amplitudes(iTE,:)/gasFit.nmrMix.amp(1);
    
    % Calculate TE 90
    deltaF(iTE) = nmrFit.nmrMix.freq(barrier_te90_idx)-nmrFit.nmrMix.freq(rbc_te90_idx);
    deltaPhase = nmrFit.nmrMix.phase(barrier_te90_idx)-nmrFit.nmrMix.phase(rbc_te90_idx);
    time180(iTE) = abs(1/(2*deltaF(iTE)));
    time90(iTE) = 0.5*time180(iTE);
    te90(iTE) = (90-deltaPhase)/(360*deltaF(iTE)) + TEs(iTE);
    relTe90 = te90(iTE) - TEs(iTE);
    while(te90(iTE)>(minte+time180(iTE)))
        % This te is too high, so subtract 180 deg of phase
        te90(iTE) = te90(iTE) - time180(iTE);
    end
    while(te90(iTE)<minte)
        % This TE is too low, so add 180 deg of phase
        te90(iTE) = te90(iTE) + time180(iTE);
    end
    
    % Show fit
    figure();
    ax = nmrFit.plotFit();
    title(ax,['TE' num2str(iTE)]);
    
    
    % Describe fit
    nmrFit.describe(gasFit.nmrMix.amp(1));
end



figure()
ax1 = subplot(4,1,1);
plot(repmat(TEs',[1 nComp]),amplitudes);
xlabel('TE');
ylabel('Signal intensity (arbs)');
legend('rbc','barrier 1','barrier 2','gas 1','gas 2')

ax2 = subplot(4,1,2);
plot(repmat(TEs',[1 nComp]), barrier_ratio);
xlabel('TE');
ylabel('Ratio with barrier');
legend('rbc:barrier','barrier 1:barrier','barrier 2:barrier',...
    'gas 1:barrier','gas 2:barrier')

ax3 = subplot(4,1,3);
plot(repmat(TEs',[1 nComp]), gas_ratio);
xlabel('TE');
ylabel('Ratio with gas');
legend('rbc:gas','barrier 1:gas','barrier 2:gas', 'gas 1:gas', 'gas 2:gas');

ax4 = subplot(4,1,4);
plot(TEs, te90*1E6);
xlabel('TE');
ylabel('TE90 (usec)');

%% Display TE90 etc
deltaF
mean_te90 = mean(te90);
stdev_te90 = std(te90);
phase90_usec = round(time90*1E6)
te90_usec = round(te90'*1E6)
barrier_ratio(:,1)'

% Sumamrize
disp(['TE90=' num2str(mean_te90*1E6) 'usec (' num2str(stdev_te90*1E6) 'usec stdev)']);
disp(['Flip angle ~' num2str(flip_angle) ' (' num2str(flip_err) ' error)']);
disp(['mean RBC:Barrier = ' num2str(mean(abs(barrier_ratio(rbc_te90_idx,:)))) ' (' num2str(std(abs(barrier_ratio(rbc_te90_idx,:)))) ' std dev)']);



