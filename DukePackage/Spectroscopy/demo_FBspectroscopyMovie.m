% Guesses
amp = [1 1 1 1];
freq = [-3860 -3840 387 -35];
fwhm = [50 35 125 200];
phase = [0 0 0 0];
zeroPad = 512;
lineBroad = 0;

% Find pfile
pfile_path = filepath('/home/scott/Desktop/subject002_065/P35840.7')

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

nToAvg = 25;
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
saveMix = [];
fits = cell(nTimePoints,1);
parfor iTimePoint = 1:nTimePoints
    startIdx = startingTimePoints(iTimePoint)
    avg_data = mean(pfile.data(:,startIdx + (1:nToAvg) - 1),2);
    
    nmrMix = NMR_Mix(avg_data,t,amp,freq,fwhm,phase,zeroPad,lineBroad,pfile.rdb.rdb_hdr_ps_mps_freq/10);
    nmrMix = nmrMix.fitTimeDomainSignal();
    nmrMix = nmrMix.sortByFreq();
    
    amp_dyn(iTimePoint,:) = nmrMix.amp(:);
    freq_dyn(iTimePoint,:) = nmrMix.freq(:);
    fwhm_dyn(iTimePoint,:) = nmrMix.fwhm(:);
    phase_dyn(iTimePoint,:) = nmrMix.phase(:);
    
    startT = t_tr(startIdx);
    stuffThis = repmat(startT + floor(0.5*nToAvg*dwell_time),[1 nComp]);
    t_dyn(iTimePoint,:) = stuffThis;
    fits{iTimePoint} = nmrMix;
end

fig = figure();

barrier_color = [    0    0.4470    0.7410];
rbc_color =     [  0.8500    0.3250    0.0980];
alveoli_color =  [ 0.4940    0.1840    0.5560];
airway_color =     [0.4660    0.6740    0.1880];

max_t = max(t_tr);
ax45lim = [0 max_t 0 1.0498e+06];
ax910lim = [0 max_t -4000 -3800];
ax1415lim = [0 max_t 0 35];
ax1920lim = [0 max_t 0 0.9]; %-0.0116

patid = '002-065';
spec_movieName = ['FB_spectro_SUBJ' patid 'spec.mp4'];
spec_vidObj = VideoWriter(spec_movieName);
spec_vidObj.FrameRate = 30;
open(spec_vidObj);

gasSpec_movieName = ['FB_spectro_SUBJ' patid 'specGas.mp4'];
gasSpec_vidObj = VideoWriter(gasSpec_movieName);
gasSpec_vidObj.FrameRate = 30;
open(gasSpec_vidObj);

disSpec_movieName = ['FB_spectro_SUBJ' patid 'specDis.mp4'];
disSpec_vidObj = VideoWriter(disSpec_movieName);
disSpec_vidObj.FrameRate = 30;
open(disSpec_vidObj);

time_movieName = ['FB_spectro_SUBJ' patid 'time.mp4'];
time_vidObj = VideoWriter(time_movieName);
time_vidObj.FrameRate = 30;
open(time_vidObj);

for iTimePoint = 1:nTimePoints
    timeStart = t_tr(startingTimePoints(iTimePoint));
    timeEnd = t_tr(startingTimePoints(iTimePoint)+nToAvg);
    disp([num2str(100*iTimePoint/nTimePoints) '% done.']);
    measured_time = fits{iTimePoint}.timeDomainSignal;
    fit_time = fits{iTimePoint}.calcTimeDomainSignal();
    measured_spec = fits{iTimePoint}.spectralDomainSignal;
    fit_spectrum = fits{iTimePoint}.calcSpectralDomainSignal();
    
    ax45 = subplot(4,5,[4 5]);
    signal_dyn = amp_dyn .* fwhm_dyn;
    plot(t_dyn(:,1),sum(signal_dyn,2),'-k','LineWidth',3);
    hold on;
    plot(t_dyn(:,1),signal_dyn(:,1),'LineWidth',2,'Color',rbc_color);
    plot(t_dyn(:,2),signal_dyn(:,2),'LineWidth',2,'Color',barrier_color);
    plot(t_dyn(:,3),signal_dyn(:,3),'LineWidth',2,'Color',airway_color);
    plot(t_dyn(:,4),signal_dyn(:,4),'LineWidth',2,'Color',alveoli_color);
    patch([timeStart timeEnd timeEnd timeStart],...
        [ax45lim(4) ax45lim(4) ax45lim(3) ax45lim(3)],'black','FaceAlpha',0.5,'EdgeColor','none')
    hold off
    ylabel('Intensity (abs)');
    xlabel('Time (sec)');
    legend('Net Signal','RBC','Barrier','Airway','Alveoli','Orientation','vertical','Location','northwest');
    axis([ax45lim]);
    
    ax910 = subplot(4,5,[9 10]);
    plot(t_dyn(:,1),freq_dyn(:,1),'LineWidth',2,'Color',rbc_color);
    hold on;
    plot(t_dyn(:,2),freq_dyn(:,2),'LineWidth',2,'Color',barrier_color);
    plot(t_dyn(:,3),freq_dyn(:,3),'LineWidth',2,'Color',airway_color);
    plot(t_dyn(:,4),freq_dyn(:,4),'LineWidth',2,'Color',alveoli_color);
    patch([timeStart timeEnd timeEnd timeStart],...
        [ax910lim(4) ax910lim(4) ax910lim(3) ax910lim(3)],'black','FaceAlpha',0.5,'EdgeColor','none')
    hold off;
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    legend('RBC','Barrier','Airway','Alveoli','Orientation','vertical','Location','northwest');
    axis([ax910lim]);
    
    ax1415 = subplot(4,5,[14 15]);
    airway_dyn = signal_dyn(:,3);
    rbc2airway = signal_dyn(:,1)./airway_dyn;
    barrier2airway = signal_dyn(:,2)./airway_dyn;
    alveoli2airway = signal_dyn(:,4)./airway_dyn;
    plot(t_dyn(:,3),rbc2airway,'LineWidth',2,'Color',rbc_color);
    hold on;
    plot(t_dyn(:,2),barrier2airway,'LineWidth',2,'Color',barrier_color);
    plot(t_dyn(:,1),alveoli2airway,'LineWidth',2,'Color',alveoli_color);
    patch([timeStart timeEnd timeEnd timeStart],...
        [ax1415lim(4) ax1415lim(4) ax1415lim(3) ax1415lim(3)],'black','FaceAlpha',0.5,'EdgeColor','none')
    hold off;
    legend('RBC:Airway','Barrier:Airway','Alveoli:Airway','Orientation','vertical','Location','northwest');
    ylabel('Airway Ratio (arbs)');
    xlabel('Time (sec)');
    axis([ax1415lim]);
    
    ax1920 = subplot(4,5,[19 20]);
    barrier_dyn = signal_dyn(:,2);
    alveoli2barrier = signal_dyn(:,4)./barrier_dyn;
    airway2barrier = signal_dyn(:,3)./barrier_dyn;
    rbc2barrier = signal_dyn(:,1)./barrier_dyn;
    plot(t_dyn(:,3),rbc2barrier,'LineWidth',2,'Color',rbc_color);
    hold on;
    plot(t_dyn(:,2),airway2barrier,'LineWidth',2,'Color',airway_color);
    plot(t_dyn(:,1),alveoli2barrier,'LineWidth',2,'Color',alveoli_color);
    patch([timeStart timeEnd timeEnd timeStart],...
        [ax1920lim(4) ax1920lim(4) ax1920lim(3) ax1920lim(3)],'black','FaceAlpha',0.5,'EdgeColor','none')
    hold off;
    legend('RBC:Barrier','Airway:Barrier','Alveoli:Barrier','Orientation','vertical','Location','northwest');
    xlabel('Time (sec)');
    ylabel('Barrier Ratio (arbs)');
    axis([ax1920lim]);
    
    %% Plot freq 
    mag_axis = [-4500 1000 0 15];
    ax1 = subplot(4,5,1);
    plot(fits{iTimePoint}.f,abs(measured_spec),'-b','LineWidth',2);
    axis(mag_axis);
    ylabel('Measured');
    xlabel('Freq (Hz)');
    title('Magnitude');
    
    ax6 = subplot(4,5,6);
    plot(fits{iTimePoint}.f,abs(sum(fit_spectrum,2)),'-k','LineWidth',2);
    axis(mag_axis);
    xlabel('Freq (Hz)');
    ylabel('Fit');
    
    ax11 = subplot(4,5,11);
    plot(fits{iTimePoint}.f,abs(measured_spec-sum(fit_spectrum,2)),'-r','LineWidth',2);
    axis(mag_axis);
    xlabel('Freq (Hz)');
    ylabel('Residual');
    
    ax16 = subplot(4,5,16);
    plot(fits{iTimePoint}.f,abs(fit_spectrum(:,1)),'LineWidth',1,'Color',rbc_color);
    hold on;
    plot(fits{iTimePoint}.f,abs(fit_spectrum(:,2)),'LineWidth',1,'Color',barrier_color);
    plot(fits{iTimePoint}.f,abs(fit_spectrum(:,3)),'LineWidth',1,'Color',airway_color);
    plot(fits{iTimePoint}.f,abs(fit_spectrum(:,4)),'LineWidth',1,'Color',alveoli_color);
    hold off;
    axis(mag_axis);
    l = legend('RBC','Barrier','Airway','Alveoli','Orientation','vertical','Location','eastoutside');
    xlabel('Freq (Hz)');
    set(l,'Position',[0.0606    0.2011    0.0528    0.0673]);
    
    % Real
    real_axis = [-4500 1000 -15 15];
    ax2 = subplot(4,5,2);
    plot(fits{iTimePoint}.f,real(measured_spec),'-b','LineWidth',2);
    axis(real_axis);
    xlabel('Freq (Hz)');
    title('Real');
    
    ax7 = subplot(4,5,7);
    plot(fits{iTimePoint}.f,real(sum(fit_spectrum,2)),'-k','LineWidth',2);
    axis(real_axis);
    xlabel('Freq (Hz)');
    
    ax12 = subplot(4,5,12);
    plot(fits{iTimePoint}.f,real(measured_spec-sum(fit_spectrum,2)),'-r','LineWidth',2);
    axis(real_axis);
    
    ax17 = subplot(4,5,17);
    plot(fits{iTimePoint}.f,real(fit_spectrum(:,1)),'LineWidth',1,'Color',rbc_color);
    hold on;
    plot(fits{iTimePoint}.f,real(fit_spectrum(:,2)),'LineWidth',1,'Color',barrier_color);
    plot(fits{iTimePoint}.f,real(fit_spectrum(:,3)),'LineWidth',1,'Color',airway_color);
    plot(fits{iTimePoint}.f,real(fit_spectrum(:,4)),'LineWidth',1,'Color',alveoli_color);
    hold off;
    axis(real_axis);
    xlabel('Freq (Hz)');
    
    % Imaginary
    imag_axis = [-4500 1000 -15 15];
    ax3 = subplot(4,5,3);
    plot(fits{iTimePoint}.f,imag(measured_spec),'-b','LineWidth',2);
    axis(imag_axis);
    title('Imaginary');
    xlabel('Freq (Hz)');
    
    ax8 = subplot(4,5,8);
    plot(fits{iTimePoint}.f,imag(sum(fit_spectrum,2)),'-k','LineWidth',2);
    axis(imag_axis);
    xlabel('Freq (Hz)');
    
    ax13 = subplot(4,5,13);
    plot(fits{iTimePoint}.f,imag(measured_spec-sum(fit_spectrum,2)),'-r','LineWidth',2);
    axis(imag_axis);
    xlabel('Freq (Hz)');
    
    ax18 = subplot(4,5,18);
    plot(fits{iTimePoint}.f,imag(fit_spectrum(:,1)),'LineWidth',1,'Color',rbc_color);
    hold on;
    plot(fits{iTimePoint}.f,imag(fit_spectrum(:,2)),'LineWidth',1,'Color',barrier_color);
    plot(fits{iTimePoint}.f,imag(fit_spectrum(:,3)),'LineWidth',1,'Color',airway_color);
    plot(fits{iTimePoint}.f,imag(fit_spectrum(:,4)),'LineWidth',1,'Color',alveoli_color);
    hold off;
    axis(imag_axis);
    xlabel('Freq (Hz)');
    
    set(ax1,'XDir','reverse');
    set(ax2,'XDir','reverse');
    set(ax3,'XDir','reverse');
    set(ax6,'XDir','reverse');
    set(ax7,'XDir','reverse');
    set(ax8,'XDir','reverse');
    set(ax11,'XDir','reverse');
    set(ax12,'XDir','reverse');
    set(ax13,'XDir','reverse');
    set(ax16,'XDir','reverse');
    set(ax17,'XDir','reverse');
    set(ax18,'XDir','reverse');
    writeVideo(spec_vidObj,getframe(fig));
    
    % Zoom in on gas peaks
    gas_xlim = [-5000 -3000];
    axes(ax1); xlim(gas_xlim);
    axes(ax2); xlim(gas_xlim);
    axes(ax3); xlim(gas_xlim);
    axes(ax6); xlim(gas_xlim);
    axes(ax7); xlim(gas_xlim);
    axes(ax8); xlim(gas_xlim);
    axes(ax11); xlim(gas_xlim);
    axes(ax12); xlim(gas_xlim);
    axes(ax13); xlim(gas_xlim);
    axes(ax16); xlim(gas_xlim);
    axes(ax17); xlim(gas_xlim);
    axes(ax18); xlim(gas_xlim);
    set(ax1,'XDir','reverse');
    set(ax2,'XDir','reverse');
    set(ax3,'XDir','reverse');
    set(ax6,'XDir','reverse');
    set(ax7,'XDir','reverse');
    set(ax8,'XDir','reverse');
    set(ax11,'XDir','reverse');
    set(ax12,'XDir','reverse');
    set(ax13,'XDir','reverse');
    set(ax16,'XDir','reverse');
    set(ax17,'XDir','reverse');
    set(ax18,'XDir','reverse');
    writeVideo(gasSpec_vidObj,getframe(fig));
    
        % Zoom in on disolved peaks
    dis_xlim = [-1500 500];
    dis_ylim_mag = [0 5];
    dis_ylim_complex = [-5 5];
    axes(ax1); xlim(dis_xlim); ylim(dis_ylim_mag);
    axes(ax2); xlim(dis_xlim);ylim(dis_ylim_complex);
    axes(ax3); xlim(dis_xlim);ylim(dis_ylim_complex);
    axes(ax6); xlim(dis_xlim); ylim(dis_ylim_mag);
    axes(ax7); xlim(dis_xlim);ylim(dis_ylim_complex);
    axes(ax8); xlim(dis_xlim);ylim(dis_ylim_complex);
    axes(ax11); xlim(dis_xlim);ylim(dis_ylim_mag);
    axes(ax12); xlim(dis_xlim);ylim(dis_ylim_complex);
    axes(ax13); xlim(dis_xlim);ylim(dis_ylim_complex);
    axes(ax16); xlim(dis_xlim);ylim(dis_ylim_mag);
    axes(ax17); xlim(dis_xlim);ylim(dis_ylim_complex);
    axes(ax18); xlim(dis_xlim);ylim(dis_ylim_complex);
    set(ax1,'XDir','reverse');
    set(ax2,'XDir','reverse');
    set(ax3,'XDir','reverse');
    set(ax6,'XDir','reverse');
    set(ax7,'XDir','reverse');
    set(ax8,'XDir','reverse');
    set(ax11,'XDir','reverse');
    set(ax12,'XDir','reverse');
    set(ax13,'XDir','reverse');
    set(ax16,'XDir','reverse');
    set(ax17,'XDir','reverse');
    set(ax18,'XDir','reverse');
    writeVideo(disSpec_vidObj,getframe(fig));
    
    
    %% Plot time
    mag_axis = [0 max(t) 0 3000];
    ax1 = subplot(4,5,1);
    plot(fits{iTimePoint}.t,abs(measured_time),'-b','LineWidth',2);
    axis(mag_axis);
    ylabel('Measured');
    xlabel('Time (sec)');
    title('Magnitude');
    
    ax6 = subplot(4,5,6);
    plot(fits{iTimePoint}.t,abs(sum(fit_time,2)),'-k','LineWidth',2);
    axis(mag_axis);
    ylabel('Fit');
    xlabel('Time (sec)');
    
    ax11 = subplot(4,5,11);
    plot(fits{iTimePoint}.t,abs(measured_time-sum(fit_time,2)),'-r','LineWidth',2);
    axis(mag_axis);
    ylabel('Residual');
    xlabel('Time (sec)');
    
    ax16 = subplot(4,5,16);
    plot(fits{iTimePoint}.t,abs(fit_time(:,1)),'LineWidth',1,'Color',rbc_color);
    hold on;
    plot(fits{iTimePoint}.t,abs(fit_time(:,2)),'LineWidth',1,'Color',barrier_color);
    plot(fits{iTimePoint}.t,abs(fit_time(:,3)),'LineWidth',1,'Color',airway_color);
    plot(fits{iTimePoint}.t,abs(fit_time(:,4)),'LineWidth',1,'Color',alveoli_color);
    hold off;
    axis(mag_axis);
    xlabel('Time (sec)');
    l = legend('RBC','Barrier','Airway','Alveoli','Orientation','vertical','Location','eastoutside');
    set(l,'Position',[0.0606    0.2011    0.0528    0.0673]);
    
    
    % Real
    real_axis = [0 max(t) -2000 3000];
    ax2 = subplot(4,5,2);
    plot(fits{iTimePoint}.t,real(measured_time),'-b','LineWidth',2);
    axis(real_axis);
    title('Real');
    xlabel('Time (sec)');
    
    ax7 = subplot(4,5,7);
    plot(fits{iTimePoint}.t,real(sum(fit_time,2)),'-k','LineWidth',2);
    axis(real_axis);
    xlabel('Time (sec)');
    
    ax12 = subplot(4,5,12);
    plot(fits{iTimePoint}.t,real(measured_time-sum(fit_time,2)),'-r','LineWidth',2);
    axis(real_axis);
    xlabel('Time (sec)');
    
    ax17 = subplot(4,5,17);
    plot(fits{iTimePoint}.t,real(fit_time(:,1)),'LineWidth',1,'Color',rbc_color);
    hold on;
    plot(fits{iTimePoint}.t,real(fit_time(:,2)),'LineWidth',1,'Color',barrier_color);
    plot(fits{iTimePoint}.t,real(fit_time(:,3)),'LineWidth',1,'Color',airway_color);
    plot(fits{iTimePoint}.t,real(fit_time(:,4)),'LineWidth',1,'Color',alveoli_color);
    hold off;
    axis(real_axis);
    xlabel('Time (sec)');
    
    % Imaginary
    imag_axis = [0 max(t) -2000 3000];
    ax3 = subplot(4,5,3);
    plot(fits{iTimePoint}.t,imag(measured_time),'-b','LineWidth',2);
    axis(imag_axis);
    title('Imaginary');
    xlabel('Time (sec)');
    
    ax8 = subplot(4,5,8);
    plot(fits{iTimePoint}.t,imag(sum(fit_time,2)),'-k','LineWidth',2);
    axis(imag_axis);
    xlabel('Time (sec)');
    
    ax13 = subplot(4,5,13);
    plot(fits{iTimePoint}.t,imag(measured_time-sum(fit_time,2)),'-r','LineWidth',2);
    axis(imag_axis);
    xlabel('Time (sec)');
    
    ax18 = subplot(4,5,18);
    plot(fits{iTimePoint}.t,imag(fit_time(:,1)),'LineWidth',1,'Color',rbc_color);
    hold on;
    plot(fits{iTimePoint}.t,imag(fit_time(:,2)),'LineWidth',1,'Color',barrier_color);
    plot(fits{iTimePoint}.t,imag(fit_time(:,3)),'LineWidth',1,'Color',airway_color);
    plot(fits{iTimePoint}.t,imag(fit_time(:,4)),'LineWidth',1,'Color',alveoli_color);
    hold off;
    axis(imag_axis);
    hold on;
    xlabel('Time (sec)');
    
    set(ax1,'XDir','normal');
    set(ax2,'XDir','normal');
    set(ax3,'XDir','normal');
    set(ax6,'XDir','normal');
    set(ax7,'XDir','normal');
    set(ax8,'XDir','normal');
    set(ax11,'XDir','normal');
    set(ax12,'XDir','normal');
    set(ax13,'XDir','normal');
    set(ax16,'XDir','normal');
    set(ax17,'XDir','normal');
    set(ax18,'XDir','normal');
    writeVideo(time_vidObj,getframe(fig));
end

% Close the movie files.
close(spec_vidObj);
close(gasSpec_vidObj);
close(disSpec_vidObj);
close(time_vidObj);




