% Guesses 
amp = [1 1 1 1 1];
freq = [-35 -285 -293 -3840 -3870];
fwhm = [215 200 115 50 50];
phase = [0 0 0 0 0];

% Find pfile
pfile_path = filepath('C:\Users\Scott\Desktop\')

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

nToAvg = 20;
skipSize = 1;

% Fit spectra
nComp = length(amp);
startingTimePoints = 1:skipSize:(nFrames-nToAvg);
nTimePoints = length(startingTimePoints);
amp_dyn = zeros(nTimePoints,nComp);
freq_dyn = zeros(nTimePoints,nComp);
fwhm_dyn = zeros(nTimePoints,nComp);
phase_dyn = zeros(nTimePoints,nComp);
t_dyn = zeros(nTimePoints,nComp);
parfor iTimePoint = 1:nTimePoints
	startIdx = startingTimePoints(iTimePoint)
	avg_data = mean(pfile.data(:,startIdx + (1:nToAvg) - 1),2);
	
	nmrMix = NMR_Mix(avg_data,t,amp,freq,fwhm,phase,[],[],pfile.rdb.rdb_hdr_ps_mps_freq/10);
%     nmrMix = NMR_Mix(avg_data,t,[],[],[],[],[],[],pfile.rdb.rdb_hdr_ps_mps_freq/10);
	nmrMix = nmrMix.fitTimeDomainSignal();
%     nmrMix = nmrMix.autoAddComponents(3);
    nmrMix = nmrMix.sortByFreq();
	
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

deltaf = [freq_dyn(:,1)-freq_dyn(300,1)...
	freq_dyn(:,2)-freq_dyn(300,2)...
	freq_dyn(:,3)-freq_dyn(300,3)...
    freq_dyn(:,4)-freq_dyn(300,4)...
    freq_dyn(:,5)-freq_dyn(300,5)];

ratios = [rbc2barrier rbc2gas barrier2gas];

figure()
ax11 = subplot(3,5,1);
plot(t_dyn(:,1),freq_dyn(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),freq_dyn(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),freq_dyn(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),freq_dyn(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),freq_dyn(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('Comp Freq (Hz)');

ax21 = subplot(3,5,6);
plot(t_dyn(:,1),deltaf(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),deltaf(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),deltaf(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),deltaf(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),deltaf(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('Diff Freq (Hz)');

ax31 = subplot(3,5,11);
plot(t_dyn(:,1),phase_dyn(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),phase_dyn(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),phase_dyn(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),phase_dyn(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),phase_dyn(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('Phase (Degrees)');


intensity =amp_dyn.*fwhm_dyn; 
ax12 = subplot(3,5,2);
plot(t_dyn(:,1),intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('Intensity(arbs)');

ax22 = subplot(3,5,7);
plot(t_dyn(:,1),amp_dyn(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),amp_dyn(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),amp_dyn(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),amp_dyn(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),amp_dyn(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('Amplitude (arbs)');

ax32 = subplot(3,5,12);
plot(t_dyn(:,1),fwhm_dyn(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),fwhm_dyn(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),fwhm_dyn(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),fwhm_dyn(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),fwhm_dyn(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('FWHM (Hz)');

ax13 = subplot(3,5,3);
multiply_factor = 1./sum(intensity,2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:Intensity');

ax23 = subplot(3,5,8);
multiply_factor = 1./sum(intensity(:,4:5),2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:Gases');

ax33 = subplot(3,5,13);
multiply_factor = 1./sum(intensity(:,1:3),2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:Dissolved');

ax14 = subplot(3,5,4);
multiply_factor = 1./sum(intensity(:,1),2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:RBC');

ax24 = subplot(3,5,9);
multiply_factor = 1./sum(intensity(:,2),2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:Fibrosis');

ax34 = subplot(3,5,14);
multiply_factor = 1./sum(intensity(:,3),2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:Barrier');

ax15 = subplot(3,5,5);
multiply_factor = 1./sum(intensity(:,4),2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:Airways');

ax25 = subplot(3,5,10);
multiply_factor = 1./sum(intensity(:,5),2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:Alveoli');

ax35 = subplot(3,5,15);
multiply_factor = 1./sum(intensity(:,2:3),2);
plot(t_dyn(:,1),multiply_factor.*intensity(:,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold on;
plot(t_dyn(:,2),multiply_factor.*intensity(:,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(:,3),multiply_factor.*intensity(:,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(t_dyn(:,4),multiply_factor.*intensity(:,4),'-','LineWidth',2,'color',[0.3010    0.7450    0.9330]);
plot(t_dyn(:,5),multiply_factor.*intensity(:,5),'-','LineWidth',2,'color',[0    0.4470    0.7410]);
hold off
xlabel('Time (sec)');
ylabel('X:(Barrier+Fibrosis)');

legend('RBC','Fibrosis','Barrier','Airways','Alveoli','Location','North','orientation','horizontal')
linkaxes([ax11 ax21 ax31 ax12 ax22 ax32 ax13 ax23 ax33 ax14 ax24 ax34 ax15 ax25 ax35],'x');


