%%3peaks
%   Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
dis_fit_guess = [
    1           -100            136          0; % Component #1
    1           -348           130          0; % Component #3
    1           -3837           35          0; % Component #4
    ];

rbc_idx = 1;
barrier_idx = 2;
rbc_te90_idx = 1;
barrier_te90_idx = 2;
gas_idx = 3;

% %%4peaks
% %   Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
% dis_fit_guess = [
%     1           -28            190          0; % Component #1
%     1           -290           200          0; % Component #2
%     1           -381           130          0; % Component #3
%     1           -3850           30          0; % Component #4
%     ];
% 
% rbc_idx = 1;
% barrier_idx = 2:3;
% rbc_te90_idx = 1;
% barrier_te90_idx = 3;
% gas_idx = 4;

% %%5 peaks
% %   Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
% dis_fit_guess = [
%     1           0            215          0; % Component #1
%     1           -280           200          0; % Component #2
%     1           -380           115          0; % Component #3
%     1           -3835           30          0; % Component #4
%     1           -3860           60          0; % Component #4
%     ];
% 
% rbc_idx = 1;
% barrier_idx = 2:3;
% rbc_te90_idx = 1;
% barrier_te90_idx = 3;
% gas_idx = 4:5;

gas_fit_guess = [   1           8           15          0]; % Component #

% TE Limits
maxTEghetto = 1500E-6;
minTEghetto = 600E-6;
minte = 400E-6;

% Miscelaneous parameters
endToff = 1;
zeropadsize = 512;
linebroadening = 0; %Hz
skipDownstreamFrames = 25;
throwAwayCalFrames = 0;
TEs = [875 975 1075 1175]*1E-6;
nDisFrames = 200;
nGasFrames = 1;
nTE = length(TEs);
nComp = size(dis_fit_guess,1);

colors = [     0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4660    0.6740    0.1880
    0    0.4470    0.7410
    0.3010    0.7450    0.9330
    0.4940    0.1840    0.5560
    0.6350    0.0780    0.1840];
linestyles = {'-','--','-.',':'};

% Get Pfile
% pfile_path = filepath('/home/scott/Desktop/')
pfile_path = filepath('/home/scott/Public/data/20150716/tardis/')

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% Pull relavent info from header
npts = pfile.rdb.rdb_hdr_frame_size;nGas = length(gas_idx);
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
gasFit.describe();
% figure()
% gasFit.plotFit();

%% Split flip calibration frames into separate pfile
flipCal_pfile = pfile;
TE90=1215.433
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
delta_rec_freq = pfile.image.user12; % freq_off in Hz

b = zeros(4,2048);
te90 = zeros(nTE,1);
if(nComp > 3)
    plasma_phase = zeros(nTE,1);
end

specFitFig = figure('Name','Spectral Fitting','Units','Pixels','Position',[66 372 820 602]);
te90Fig = figure('Name','TE 90','Units','Pixels','Position',[887 372 494 600]);
axToLink = [];
legendstr = {'RBC','Barrier 1','Barrier 2','Gas 1' 'Gas 2'};
relTimeVals = (linspace(-6000,max(t*1E6),50000))*1E-6;
nDis = length(rbc_idx)+ length(barrier_idx);
nGas = length(gas_idx);
legStr = cell(1,(nDis+1));
metricOfGoodness = zeros(0,length(relTimeVals));
for iLeg = 1:nDis
    legStr{iLeg} = [legendstr{iLeg} ] ;
end
legStr{nDis+1} = ['Net Signal'];
legendStrings = cell(1, nDis+1);
for iComp=1:nDis
    legendStrings{iComp} = legendstr{iComp};
end
for iGas=1:nGas
    legendStrings{nDis+iGas} = legendstr{3+iGas};
end
for iTE = 1:nTE
    % Calculate average fid for TE
    teData = dis_pfile.data(:,(skipDownstreamFrames*nTE+iTE):nTE:end);
    teData = mean(teData,2);
    
    % Force zero phase at latest TE to align all FIDs
    nmrFit = NMR_Fit(teData, t, zeropadsize,linebroadening, dis_fit_guess(:,1),dis_fit_guess(:,2),...
        dis_fit_guess(:,3),dis_fit_guess(:,4));
    nmrFit = nmrFit.fitTimeDomainSignal();  
%     figure();
%     nmrFit.plotFit();
%     figure()
%     nmrFit.plotTimeFit();
    startingVec = nmrFit.nmrMix.calcTimeDomainSignal(TEs(nTE)-TEs(iTE));
    startingPhase = angle(startingVec);
    teData = teData.*exp(-1i*startingPhase);
    nmrFit = NMR_Fit(teData, t, zeropadsize,linebroadening, nmrFit.nmrMix.amp,nmrFit.nmrMix.freq,...
        nmrFit.nmrMix.fwhm,nmrFit.nmrMix.phase-rad2deg(startingPhase));
    figure();
    nmrFit.plotFit();
    
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
    
    % Calculate TE_ghetto which always has 
    onRbcResMix = NMR_Mix(nmrFit.nmrMix.amp([rbc_idx barrier_idx gas_idx]), ...
        nmrFit.nmrMix.freq([rbc_idx barrier_idx gas_idx])-nmrFit.nmrMix.freq(rbc_idx), ... % differential freq
        nmrFit.nmrMix.fwhm([rbc_idx barrier_idx gas_idx]), ...
        nmrFit.nmrMix.phase([rbc_idx barrier_idx gas_idx])-nmrFit.nmrMix.phase(rbc_idx) ... % zero starting phase 
        );
    
    relTimeTE = relTimeVals -  TEs(iTE);
    onRbcSig = onRbcResMix.calcTimeDomainSignal(relTimeTE(:));
    onRbcCompSig = onRbcResMix.calcComponentTimeDomainSignal(relTimeTE(:));
    metricOfGoodness(iTE,:) = real(onRbcCompSig(:,rbc_idx)./sum(abs(real(onRbcCompSig(:,[barrier_idx gas_idx]))),2));
    
    actSig = nmrFit.nmrMix.calcTimeDomainSignal(relTimeTE(:));
    actCompSig = nmrFit.nmrMix.calcComponentTimeDomainSignal(relTimeTE(:));
    
    % Show Time signal
    figure(te90Fig);
    ax1 = subplot(3,1,1);
    hold on;
    line((relTimeVals)*1E6,real(actCompSig(:,rbc_idx)),'Color',colors(1,:),'Linestyle',linestyles{iTE},'Parent',ax1);
    for ibar = 1:length(barrier_idx)
        line((relTimeVals)*1E6,real(actCompSig(:,1+ibar)),'Color',colors(1+ibar,:),'Linestyle',linestyles{iTE},'Parent',ax1);
    end
    line((relTimeVals)*1E6,real(actSig),'Color','k','Linestyle',linestyles{iTE},'Parent',ax1);
    hold off;
    xlabel('Time (usec)');
    ylabel('Real Signal');
    if(iTE==nTE)
        legend(legStr,'Orientation','Vertical');
    end

    % Show phase relative to RBC
    ax2 = subplot(3,1,2);
    hold on;
    line((relTimeVals)*1E6,rad2deg(angle(onRbcCompSig(:,rbc_idx))),'Color',colors(1,:),'Linestyle',linestyles{iTE},'Parent',ax2);
    for ibar = 1:length(barrier_idx)
        line((relTimeVals)*1E6,rad2deg(angle(onRbcCompSig(:,1+ibar))),'Color',colors(1+ibar,:),'Linestyle',linestyles{iTE},'Parent',ax2);
    end
    line((relTimeVals)*1E6,rad2deg(angle(onRbcSig)),'Color','k','Linestyle',linestyles{iTE},'Parent',ax2);
    hold off;
    xlabel('Time (usec)');
    ylabel('Relative Phase (Degrees)');
    if(iTE==nTE)
        legend(legStr,'Orientation','Vertical');
    end
    ylim(ax2,[-180 180]);
    set(ax2,'YTick',180*[-1 -0.5 0 0.5 1],'YGrid','on','GridLineStyle',':');
    
    % Show metricOfGoodness
    ax3 = subplot(3,1,3);
    hold on;
    line((relTimeVals)*1E6,metricOfGoodness(iTE,:),'Linestyle',linestyles{iTE},'Parent',ax3);
    hold off;
    xlabel('Time (usec)')
    ylabel('RBC:barrierContamination ratio')
    linkaxes([ax1; ax2; ax3],'x');
    xlim(ax3,[0 3000]);
    ylim(ax3,[0 5]);
    
    % Show component fits
    figure(specFitFig);
    ax3 = subplot(nTE,2,1+2*((1:nTE)-1));
    hold on;
    if(iTE == 1)
        axToLink = [axToLink ax3];
    end
    nSamples = (2^16-1);
    dwell_time = nmrFit.t(2)-nmrFit.t(1);
    fineFreq = linspace(-0.5,0.5,nSamples)/dwell_time;
    nComponents = length(nmrFit.nmrMix.amp);
    scaledAmp = nmrFit.nmrMix.amp;
    deltaT = TEs(iTE)-TEs(1);
    for iComp = 1:nComponents
        scaledAmp(iComp) = nmrFit.nmrMix.amp(iComp)*exp(deltaT*pi*nmrFit.nmrMix.fwhm(iComp));
    end
    scaledNmrMix = NMR_Mix(scaledAmp,nmrFit.nmrMix.freq,...
        nmrFit.nmrMix.fwhm,nmrFit.nmrMix.phase);
    
    
    individualSpectrums = scaledNmrMix.calcComponentLorentzianCurves(fineFreq(:));
    
    
    for iComp=1:nDis
        line(fineFreq,real(individualSpectrums(:,iComp)),'Color',colors(iComp,:),'Linestyle',linestyles{iTE},'Parent',ax3);
    end
    for iGas=1:nGas
       line(fineFreq,real(individualSpectrums(:,nDis+iGas)),'Color',colors(3+iGas,:),'Linestyle',linestyles{iTE},'Parent',ax3);
    end
    
    if(iTE==nTE)
        legend(legendStrings,'Location','North','Orientation','vertical');
    end 
    xlabel('Freq(Hz)');
    ylabel('Lorentzian Intensity');
    set(ax3,'XDir','reverse');
    hold off;
    
    % Show residuals at first TE
    ax4 = subplot(nTE,2,2*(iTE-1)+2);
    axToLink = [axToLink ax4];  
    fittedSpectrum = nmrFit.nmrMix.calcSpectralDomainSignal(nmrFit.f);
    residualSpectrum = nmrFit.spectralDomainSignal - fittedSpectrum;
    axToLink = [axToLink ax4];
    line(nmrFit.f,real(nmrFit.spectralDomainSignal),'Color','b','Linestyle','-','Parent',ax4);
    hold on;
    line(nmrFit.f,real(fittedSpectrum),'Color','g','Linestyle','-','Parent',ax4);
    line(nmrFit.f,real(residualSpectrum),'Color','r','Linestyle','-','Parent',ax4);
    hold off;
    if(iTE==1)
        legend('Measured','Fitted','Residual','Location','North','Orientation','horizontal');
    end
    xlabel('Freq (Hz)');
    ylabel('Real Intensity');
    set(ax4,'XDir','reverse');
    title(ax4,['TE' num2str(iTE)]);
    linkaxes(axToLink,'x');
    xlim(ax4,[-725 250]);
    
    % Describe fit
    nmrFit.describe(gasFit.nmrMix.amp(1));
end
worstCaseMetric = min(metricOfGoodness,[],1);
achievableTime = (relTimeVals > minTEghetto) & (relTimeVals < maxTEghetto);
worstAchievable = worstCaseMetric(achievableTime);
timeAchievable = relTimeVals(achievableTime);
[rbc_to_crap optimum_idx] = max(worstAchievable);
te_ghetto = timeAchievable(optimum_idx);
figure(te90Fig);
ax3 = subplot(3,1,3);
hold on;
plot(te_ghetto*1E6,rbc_to_crap,'o','MarkerEdgeColor','none','MarkerFaceColor','r');
hold off;

figure('name','Spectroscopic Ratios','Units','Pixels','Position',[1383 371 538 603])
ax1 = subplot(4,1,1);
plot(TEs',amplitudes(:,1),'Color',colors(1,:));
hold on;
for iDis = 2:nDis
    plot(TEs',amplitudes(:,iDis),'Color',colors(iDis,:));
end
for iGas = 1:nGas
    plot(TEs',amplitudes(:,nDis+iGas),'Color',colors(3+nGas,:));
end
hold off;
xlabel('TE');
ylabel('Signal intensity (arbs)');
legend(legendStrings);

ax2 = subplot(4,1,2);
bar_leg = {'rbc:barrier','barrier 1:barrier','barrier 2:barrier',...
    'gas 1:barrier','gas 2:barrier'};
plot(TEs',barrier_ratio(:,1),'Color',colors(1,:));
hold on;
barrierlegendStrings = cell(1, nDis+nGas);
barrierlegendStrings{1} = bar_leg{1};
for iDis=2:nDis
    plot(TEs',barrier_ratio(:,iDis),'Color',colors(iDis,:));
    barrierlegendStrings{iDis} = bar_leg{iDis};
end
for iGas=1:nGas
    plot(TEs',barrier_ratio(:,nDis+iGas),'Color',colors(3+nGas,:));
    barrierlegendStrings{nDis+iGas} = bar_leg{3+iGas};
end
hold off;
xlabel('TE');
ylabel('Ratio with barrier');
legend(barrierlegendStrings)

ax3 = subplot(4,1,3);
gas_leg = {'rbc:gas','barrier 1:gas','barrier 2:gas', 'gas 1:gas', 'gas 2:gas'};
plot(TEs',gas_ratio(:,1),'Color',colors(1,:));
hold on;
gasLegendStrings = cell(1, nDis+nGas);
gasLegendStrings{1} = gas_leg{1};
for iDis=2:nDis
    plot(TEs',gas_ratio(:,iDis),'Color',colors(iDis,:));
    gasLegendStrings{iDis} = gas_leg{iDis};
end
for iGas=1:nGas
    plot(TEs',gas_ratio(:,nDis+iGas),'Color',colors(3+nGas,:));
    gasLegendStrings{nDis+iGas} = gas_leg{3+iGas};
end
hold off;
xlabel('TE');
ylabel('Ratio with gas');
legend(gasLegendStrings)


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
rbcToBarrierRatio = barrier_ratio(:,1)'

% Sumamrize
disp(['TE_ghetto' '=' num2str(te_ghetto*1E6)]);
disp(['TE90=' num2str(mean_te90*1E6) 'usec (' num2str(stdev_te90*1E6) 'usec stdev)']);
disp(['Flip angle ~' num2str(flip_angle) ' (' num2str(flip_err) ' error)']);
disp(['mean RBC:Barrier = ' num2str(mean(abs(barrier_ratio(:,rbc_te90_idx)))) ' (' num2str(std(abs(barrier_ratio(:,rbc_te90_idx)))) ' std dev)']);

te_ghetto = round(te_ghetto*1E6)
rbc_to_crap

