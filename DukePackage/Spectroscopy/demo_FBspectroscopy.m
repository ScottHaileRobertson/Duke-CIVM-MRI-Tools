nToAvg = 15;
skipSize = 1;
startInhale = 1;

rbc_idx = 1;
barrier_idx = 2:3;
gas_idx = [4];

% Fit middle third of data
% area = [1 1 1 1 1];
% freq = [-35 -285 -393 -3840 -3870];
% fwhm = [215 200 115 50 50];
% phase = [0 0 0 0 0];
area_orig = [1 1 1 1 1];
freq_orig = [-35 -285 -393 -3765 -3842 ];
fwhm_orig = [215 200 150 70 30];
phase_orig = [0 0 0 0 0 ];

% Find pfile
pfile_path = filepath('/home/scott/Desktop/Subj65/P35840_fbspect.7')

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

avgdata = mean(pfile.data(:,300:(end-201)),2);
nmrFit = NMR_TimeFit(avgdata, t, area_orig,freq_orig,fwhm_orig,phase_orig,[],[]);
nmrFit = nmrFit.fitTimeDomainSignal();

% Create guesses from fitted data
area = nmrFit.area;
freq = nmrFit.freq;
fwhm = nmrFit.fwhm;
phase = nmrFit.phase;

% Upper bounds
freq_halfRange = [25 25 25 2 20 ]; 
fwhm_halfRange = [30 100 100 30 15];
area_ub = inf*ones(size(fwhm));
freq_ub = nmrFit.freq + freq_halfRange;
fwhm_ub = nmrFit.fwhm + fwhm_halfRange;
phase_ub = inf*ones(size(fwhm));

% Lower Bounds
area_lb = zeros(size(area));
freq_lb = nmrFit.freq - freq_halfRange;
fwhm_lb = max(0,nmrFit.fwhm - fwhm_halfRange);
phase_lb = -inf*ones(size(fwhm));

% % Override gas peak fits
% freq_ub(4:5) = freq_orig(4:5) + freq_halfRange(4:5);
% fwhm_ub(4:5) = fwhm_orig(4:5) + fwhm_halfRange(4:5);
% freq_lb(4:5) = freq_orig(4:5) - freq_halfRange(4:5);
% fwhm_lb(4:5) = max(0,fwhm_orig(4:5) - fwhm_halfRange(4:5));
% 
% % freq = freq_orig
% freq = [6 -295 -365 -3765 -3840];
% fwhm = [195 200 154 80 19.2];
% phase = [80 -58 18 -25 -70];
% 
% % Override all fits/guesses
% % freq_lb = [6 -295 -365 -3765 -3840];
% % fwhm_lb = [160 50 50 50 10];
% % phase_lb = [60 -200 -10 -80 -85];
% % 
% % freq_ub = [50 -250 -355 -3760 -3837];
% % fwhm_ub = [240 500 250 80 30];
% % phase_ub = [100 0 50 50 -50];
% freq_lb = freq;
% fwhm_lb = fwhm;
% phase_lb = phase;
% 
% freq_ub = freq;
% fwhm_ub = fwhm;
% phase_ub = phase;

% Fit spectra
nComp = length(area);
startingTimePoints = 1:skipSize:(nFrames-nToAvg);
nTimePoints = length(startingTimePoints);
area_dyn = zeros(nTimePoints,nComp);
freq_dyn = zeros(nTimePoints,nComp);
fwhm_dyn = zeros(nTimePoints,nComp);
phase_dyn = zeros(nTimePoints,nComp);
t_dyn = zeros(nTimePoints,nComp);
tic
parfor iTimePoint = 1:nTimePoints
    startIdx = startingTimePoints(iTimePoint)
	avgdata = mean(pfile.data(:,startIdx + (1:nToAvg) - 1),2);     
        
      nmrFit = NMR_TimeFit(avgdata, t, area,freq,fwhm,phase,[],[]);
%     nmrFit = nmrFit.setBounds( area_lb, area_ub, freq_lb, freq_ub,...
%         fwhm_lb, fwhm_ub, phase_lb, phase_ub);
    nmrFit = nmrFit.fitTimeDomainSignal();
%     nmrFit.plotFit()
    
    
    area_dyn(iTimePoint,:) = nmrFit.area(:);
    freq_dyn(iTimePoint,:) = nmrFit.freq(:);
    fwhm_dyn(iTimePoint,:) = nmrFit.fwhm(:);
    phase_dyn(iTimePoint,:) = nmrFit.phase(:);
    
    startT = t_tr(startIdx);
    stuffThis = repmat(startT + floor(0.5*nToAvg*dwell_time),[1 nComp]);
    t_dyn(iTimePoint,:) = stuffThis;

    
end
toc

iMiddle = round(0.5*length(freq_dyn));
deltaf = [freq_dyn(:,1)-freq_dyn(iMiddle,1)...
    freq_dyn(:,2)-freq_dyn(iMiddle,2)...
    freq_dyn(:,3)-freq_dyn(iMiddle,3)...
    freq_dyn(:,4)-freq_dyn(iMiddle,4)];

% Colors
colors = [     0.8500    0.3250    0.0980 
    0.9290    0.6940    0.1250
    0.4660    0.6740    0.1880
     0    0.4470    0.7410
     0.3010    0.7450    0.9330
     0.4940    0.1840    0.5560
     0.6350    0.0780    0.1840];

 linewidth = 1;
titles = {'RBC' 'Barrier 1' 'Barrier 2' 'Gas 1' 'Gas 2'};
figure();
axHandles = zeros(1,4*nComp+4);
t_plot = t_dyn(startInhale:end,1)-t_dyn(startInhale,1);
for iComp = 1:nComp
    % Area
    axHandles(nComp*(iComp-1)+1) = subplot(4,nComp,iComp);
    plot(t_plot,area_dyn(startInhale:end,iComp)/gasFit.area(1),'-','LineWidth',linewidth,'color',colors(iComp,:));
    xlabel('Time (sec)');
    ylabel('Intensity (rel to ded. Gas)');
    title(titles(iComp));
    
    % Frequency
    axHandles(nComp*(iComp-1)+2) = subplot(4,nComp,nComp+iComp);
    plot(t_plot,freq_dyn(startInhale:end,iComp),'-','LineWidth',linewidth,'color',colors(iComp,:));
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    
    % Linewidth
    axHandles(nComp*(iComp-1)+3) = subplot(4,nComp,2*nComp+iComp);
    plot(t_plot,fwhm_dyn(startInhale:end,iComp),'-','LineWidth',linewidth,'color',colors(iComp,:));
    xlabel('Time (sec)');
    ylabel('Linewidth (Hz)');
    
    % Phase
    axHandles(nComp*(iComp-1)+4) = subplot(4,nComp,3*nComp + iComp);
    plot(t_plot,phase_dyn(startInhale:end,iComp),'-','LineWidth',linewidth,'color',colors(iComp,:));
    xlabel('Time (sec)');
    ylabel('Phase (degrees)');
end
linkaxes(axHandles,'x');

figure()
% Dissolved Ratio
axHandles(4*nComp + 1) = subplot(4,1,1);
multiply_factor = 1./sum(area_dyn(startInhale:end,[rbc_idx barrier_idx]),2);
plot(t_dyn(startInhale:end,3)-t_dyn(startInhale,1),area_dyn(startInhale:end,1)./area_dyn(startInhale:end,3),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
xlabel('Time (sec)');
ylabel('RBC:Barrier');

axHandles(4*nComp + 2) = subplot(4,1,2);
multiply_factor = 1./sum(area_dyn(startInhale:end,[rbc_idx barrier_idx]),2);
plot(t_dyn(startInhale:end,3)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
hold on;
plot(t_dyn(startInhale:end,2)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(startInhale:end,1)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold off
xlabel('Time (sec)');
ylabel('Ratio with Dissolved');

% Gas Ratio
axHandles(4*nComp + 3) = subplot(4,1,3);
multiply_factor = 1./sum(area_dyn(startInhale:end,[gas_idx]),2);
plot(t_dyn(startInhale:end,3)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
hold on;
plot(t_dyn(startInhale:end,2)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(startInhale:end,1)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold off
xlabel('Time (sec)');
ylabel('Ratio with Gas');

% Dedicated Gas Ratio
axHandles(4*nComp + 4) = subplot(4,1,4);
multiply_factor = 1./gasFit.area(1);
plot(t_dyn(startInhale:end,3)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,3),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
hold on;
plot(t_dyn(startInhale:end,2)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,2),'-','LineWidth',2,'color',[0.9290    0.6940    0.1250]);
plot(t_dyn(startInhale:end,1)-t_dyn(startInhale,1),multiply_factor.*area_dyn(startInhale:end,1),'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
hold off
xlabel('Time (sec)');
ylabel('Ratio with Dedicated Gas');

% Link all axes
linkaxes(axHandles,'x');
xlim([0 max(t_plot(:))]);
legend('Barrier 1','Barrier 1','RBC');


