function meanRbc2barrier = demo_dixon_DynamicSpectAsPhaseCal(varargin)

if(nargin < 1 | ~exist(varargin{1}))
    disp('Select Phase Calibration pfile');
    pfile_path = filepath('/home/scott/Public/data/')
else
    pfile_path = varargin{1};
end

nToAvg = 15;
skipSize = 1;
startInhale = 1;

rbc_idx = 1;
barrier_idx = 2:3;
gas_idx = [4];

% Fit guesses
area_orig = [1 1 1];
freq_orig = [-35 -393 -3765 ];
fwhm_orig = [215 150 70];
phase_orig = [0 0 0];

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% Create pfile for gas frame
gas_pfile = pfile;
gas_pfile.data = pfile.data(:,end);
gas_pfile.rdb.rdb_hdr_user20 = 1;

% Create pfile for dissolved frames
pfile.data = pfile.data(:,1:(end-1));
pfile.rdb.rdb_hdr_user20 = pfile.rdb.rdb_hdr_user20 -1;

% Create array of sample times (sec)
npts = pfile.rdb.rdb_hdr_frame_size;                   % Number of samples
bw = 1000*pfile.rdb.rdb_hdr_user12;                         % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw);                                 % Time between each sample
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
t = dwell_time*(0:(npts-1))';
tr = 20/1000;%pfile.image.tr;
nFrames = pfile.rdb.rdb_hdr_user20;
t_tr = tr*((1:nFrames)-1);

% Fit gas sample
gasFit = NMR_TimeFit(gas_pfile.data, t,1,-3832,20,0, [],[]);
gasFit = gasFit.fitTimeDomainSignal();

% Just take some center frames and fit them
avgdata = mean(pfile.data(:,300:450),2);
nmrFit = NMR_TimeFit(avgdata, t, area_orig,freq_orig,fwhm_orig,phase_orig,[],[]);
nmrFit = nmrFit.fitTimeDomainSignal();

minte = 700;

% Calculate RBC:barrier, TE90, etc
rbc_te90_idx = 1;
barrier_te90_idx = 2;
meanRbc2barrier = nmrFit.area(1)/nmrFit.area(2);
deltaF = nmrFit.freq(barrier_te90_idx)-nmrFit.freq(rbc_te90_idx);
deltaPhase = nmrFit.phase(barrier_te90_idx)-nmrFit.phase(rbc_te90_idx);
time180 = abs(1/(2*deltaF));
time90 = 0.5*time180;
te90 = (90-deltaPhase)/(360*deltaF) + pfile.image.te;
while(te90<minte)
    % This TE is too low, so add 180 deg of phase
    te90 = te90 + time180;
end
while(te90>(minte+time180))
    % This te is too high, so subtract 180 deg of phase
    te90 = te90 - time180;
end

disp(['TE 90 =' num2str(te90)]);
disp(['RBC:barrier =' num2str(meanRbc2barrier)]);
