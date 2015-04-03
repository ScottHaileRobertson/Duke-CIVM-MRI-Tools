% Guesses 
amp = [0.1 1 1];
freq = [3800 371 0];
fwhm = [20 125 50];
phase = [0 0 0];

% Find pfile
pfile_path = filepath('/home/scott/Public/data/20150311/subject_002_064/P03584.7')

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% Create array of sample times (sec)
npts = pfile.rdb.rdb_hdr_frame_size;                   % Number of samples
bw = 1000*pfile.rdb.rdb_hdr_user12;                         % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw);                                 % Time between each sample
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
t = dwell_time*(0:(npts-1));
tr = 20/1000;%pfile.image.tr;
nFrames = pfile.rdb.rdb_hdr_user20;
t_tr = tr*((1:nFrames)-1);

nToAvg = 100;
skipSize = 50;

% Fit spectra
nComp = length(amp);
startingTimePoints = 1:skipSize:(nFrames-nToAvg);
nTimePoints = length(startingTimePoints);
amp_dyn = zeros(nTimePoints,nComp);
freq_dyn = zeros(nTimePoints,nComp);
fwhm_dyn = zeros(nTimePoints,nComp);
phase_dyn = zeros(nTimePoints,nComp);
t_dyn = zeros(nTimePoints,nComp);
for iTimePoint = 1:nTimePoints
	startIdx = startingTimePoints(iTimePoint)
	avg_data = mean(pfile.data(:,startIdx + (1:nToAvg) - 1),2);
	
	nmrMix = NMR_Mix(amp,freq,fwhm,phase,pfile.rdb.rdb_hdr_ps_mps_freq/10);
	nmrMix = nmrMix.fitTimeDomainSignal(avg_data,t);
	
	amp_dyn(iTimePoint,:) = nmrMix.amp(:);
	freq_dyn(iTimePoint,:) = nmrMix.freq(:);
	fwhm_dyn(iTimePoint,:) = nmrMix.fwhm(:);
	phase_dyn(iTimePoint,:) = nmrMix.phase(:);
	
	startT = t_tr(startIdx);
	stuffThis = repmat(startT + floor(0.5*nToAvg*dwell_time),[1 nComp]);
	t_dyn(iTimePoint,:) = stuffThis;
end

rbc2barrier = amp_dyn(:,3)./amp_dyn(:,2);
rbc2gas = amp_dyn(:,3)./amp_dyn(:,1);
barrier2gas = amp_dyn(:,2)./amp_dyn(:,1);

deltaf = [freq_dyn(:,1)-freq_dyn(1,1)...
	freq_dyn(:,2)-freq_dyn(1,2)...
	freq_dyn(:,3)-freq_dyn(1,3)];

ratios = [rbc2barrier rbc2gas barrier2gas];

figure()
subplot(5,1,1);
plot(t_dyn,ratios);
xlabel('Time (sec)');
ylabel('Component Amplitude (arbs)');
legend('RBC/Barrier','RBC/Gas','Barrier/Gas');

subplot(5,1,2);
plot(t_dyn,amp_dyn);
xlabel('Time (sec)');
ylabel('Component Amplitude (arbs)');

subplot(5,1,3);
plot(t_dyn,deltaf);
xlabel('Time (sec)');
ylabel('Component Frequency (Hz)');

ax4 = subplot(5,1,4);
plot(t_dyn,fwhm_dyn);
xlabel('Time (sec)');
ylabel('Component FWHM (Hz)');

a4p = get(ax4,'Position');

ax5=subplot(5,1,5);
plot(t_dyn,phase_dyn);
xlabel('Time (sec)');
ylabel('Component Phase (degrees)');
legend('Gas','Barrier','RBC','Location','SouthOutside','orientation','horizontal')
a5p = get(ax5,'Position');
set(ax5, 'Position', [a5p(1) 0.14 a5p(3) 0.108])
