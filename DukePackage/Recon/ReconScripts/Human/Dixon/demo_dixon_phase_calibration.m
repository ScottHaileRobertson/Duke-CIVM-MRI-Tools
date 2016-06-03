function rbcBarrier_3peak_te90 = demo_dixon_phase_calibration(varargin)
% There are 2 uses for this script:
%    1. (During Scan) Calibrate flip angle, TE90
%    2. Calculate RBC:barrier ratio at Dixon TE (should be near TE90)
if(nargin < 1 | ~exist(varargin{1}))
    disp('Select Phase Calibration pfile');
    phaseCal_pfile = filepath('D:\Papers\Spectroscopy_3humanPeaks\Pfiles\Test\P03584.7');
    %     phaseCal_pfile = filepath('D:\Papers\Spectroscopy_3humanPeaks\Pfiles\Healthy\P_002-0039B_NEW_CAL_05632.7')
    te_rbcBarrier = [];
else
    phaseCal_pfile = varargin{1};
    if(nargin < 2)
        te_rbcBarrier = [];
    else
        te_rbcBarrier = varargin{2}; % Given in units of seconds
    end
end

% Miscelaneous parameters
constrainRBCpeak = 1;
minte = 400E-6;
endToff = 1;
skipDownstreamFrames = 25;
throwAwayCalFrames = 0;
TEs = [875 975 1075 1175]*1E-6;
nDisFrames = 200;
nGasFrames = 1;
nTE = length(TEs);

% Constants
xeGam = 17.660445;
ox_shift = (4.5/5.5)*(0.13)*(273/298)*0.917;
xe_shift = (0.3/5.5)*(273/298)*(0.548);
susc_shift = -9.06/3;
tot_shift = -xe_shift;
tot_shift_hz = tot_shift*xeGam;

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(phaseCal_pfile);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% Pull relavent info from header
npts = pfile.rdb.rdb_hdr_frame_size;
bw = pfile.rdb.rdb_hdr_user12;                 % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw*1000);
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
pfile.rdb.rdb_hdr_user12 = 1/(2*dwell_time*1000);
t = dwell_time*(0:(npts-1)); %sec

% LINEBROADING
linebroadening = 0; %Hz
zeropadsize = 10000; % Just for display, does not affect fitting

% Remove toff points
t = t(endToff:end)-t(endToff);
pfile.data = pfile.data(endToff:end,:);
pfile.rdb.rdb_hdr_frame_size = size(pfile.data,1);
npts = pfile.rdb.rdb_hdr_frame_size;

%% Split flip calibration frames into separate pfile
flipCal_pfile = pfile;
flipCal_pfile.data = pfile.data(:,(nDisFrames+nGasFrames+throwAwayCalFrames+1):end);
nFlipCal = size(flipCal_pfile.data,2);
flipCal_pfile.rdb.rdb_hdr_user20 = nFlipCal; % nframes

% Pick a reasonable DC sample
dc_sample_idx = 5;

% Calculate flip angle
[flip_angle, flip_err] = MRI.DataProcessing.calcFlipAngle(flipCal_pfile, dc_sample_idx);

fprintf (['\nFlip angle ~' num2str(flip_angle) ' (' num2str(flip_err) ' error)\n']);

%% Fit gas-phase data to 2 peaks (airway, alveolar) for a proper reference frequency
gas_fit_guess = [
    1           0           36.0          0; % Component #1
    1           -10           49.7          0; % Component #2
    ];
gas_pfile_data = pfile.data(:,nDisFrames + (1:nGasFrames));

% Fit dedicated gas
gasFit_hz= NMR_TimeFit(gas_pfile_data,t+TEs(1),...
    gas_fit_guess(:,1),gas_fit_guess(:,2),...
    gas_fit_guess(:,3),gas_fit_guess(:,4),...
    linebroadening,zeropadsize);
gasFit_hz= gasFit_hz.fitTimeDomainSignal();

% Make two peaks reasonable in frequency
badGasFreq_low = gasFit_hz.freq < -200;
badGasFreq_high =  gasFit_hz.freq > 100;
freqChange = 70;
while(any(badGasFreq_low) || any(badGasFreq_high))
    if(any(badGasFreq_low))
        if(sum(badGasFreq_low) ==2)
            error('Coulndt find freq! All too low...');
        else
            if(badGasFreq_low(2))
                gas_fit_guess(1,2) = gasFit_hz.freq(1);
                gas_fit_guess(2,2) = gasFit_hz.freq(1)-freqChange;
            else
                error('Impossible! to have higher freq too low')
            end
        end
    else
        if(sum(badGasFreq_high)==2)
            error('Couldnt find freqs! all too high...')
        else
            if(badGasFreq_high(2))
                error('Impossible! to have lower freq too high')
            else
                gas_fit_guess(1,2) = gasFit_hz.freq(2)+freqChange;
                gas_fit_guess(2,2) = gasFit_hz.freq(2);
            end
        end
    end
    
    % Refit
    gasFit_hz= NMR_TimeFit(gas_pfile_data,t+TEs(1),...
        gas_fit_guess(:,1),gas_fit_guess(:,2),...
        gas_fit_guess(:,3),gas_fit_guess(:,4),...
        linebroadening,zeropadsize);
    gasFit_hz= gasFit_hz.fitTimeDomainSignal();
    
    freqChange = freqChange-5;
    if(freqChange <= 0)
        error('Crap, cant find 2 gas-phase peaks!');
    end
    
    % Look for bad gas freqs
    badGasFreq_low = gasFit_hz.freq < -200;
    badGasFreq_high =  gasFit_hz.freq > 100;
end

% Put frequencies into ppm units
ref_freq_gas = gasFit_hz.freq(1) + tot_shift_hz;
ref_freq = gasFit_hz.freq(1) - 3832*1.67*0.59604648 + tot_shift_hz;
gasFit_hz.freq = (gasFit_hz.freq - ref_freq_gas)/xeGam;
gasFit_hz.fwhm = gasFit_hz.fwhm/xeGam;

% Report gas-phase fits
fprintf ('\n*** 3-Peak fit of gas-phase signal ***\n');
displayPPMFits(gasFit_hz);

%% Decompose multi-echo disolved-phase data
dis_pfile = pfile;
dis_pfile_data = pfile.data(:,1:nDisFrames);

% Consider each TE separately
goodFit = zeros(nTE,1);
teFits_3peak = cell(nTE,1);
teFits_5peak = cell(nTE,1);
for iTE = 1:nTE
    % Calculate average fid for TE
    teData = dis_pfile_data(:,(skipDownstreamFrames*nTE+iTE):nTE:end);
    teData = mean(teData,2);
    
    % Adapt guesses to gas-phase reference
    dis_fit_5guess = [
        1           ref_freq+216.6*xeGam            8.4*xeGam          0; % RBC
        1           ref_freq+201.6*xeGam           15.0*xeGam          0; % B1
        1           ref_freq+196.5*xeGam           9.3*xeGam          0; % B2
        1           ref_freq+0.03*xeGam          1.9*xeGam          0; % Airways
        1           ref_freq-2.3*xeGam          3.7*xeGam          0; %Alveoli
        ];
    
    % Perform 5-peak fit
    dis5peakFit = NMR_TimeFit(teData, t+TEs(iTE), ...
        dis_fit_5guess(:,1),dis_fit_5guess(:,2),...
        dis_fit_5guess(:,3),dis_fit_5guess(:,4),...
        linebroadening, zeropadsize);
    dis5peakFit = dis5peakFit.fitTimeDomainSignal();
    
    % Store fit
    teFits_5peak{iTE} = dis5peakFit;
    
    
    % Create guesses for 3 peak fit
    dis_fit_3guess = [
        1           ref_freq+216.2*xeGam            9.0*xeGam          0; % RBC
        1           ref_freq+197.5*xeGam           8.0*xeGam          0; % Barrier
        1           ref_freq-0.2*xeGam          2.0*xeGam          0; % GAs
        ];
    if(constrainRBCpeak)
        % Perform 3-peak fit after constraining RBC freq
        dis_fit_3guess(1,2) = dis5peakFit.freq(1);
        rbcArea = dis5peakFit.area(1);
        barrierArea = sum(dis5peakFit.area(2:3));
        gasArea = sum(dis5peakFit.area(4:5));
        dis_fit_3guess(1,1) = rbcArea;
        dis_fit_3guess(2,1) = barrierArea;
        dis_fit_3guess(3,1) = gasArea;
    end
    
    % Fit to 3 peaks
    dis3peakFit = NMR_TimeFit(teData, t+TEs(iTE), ...
        dis_fit_3guess(:,1),dis_fit_3guess(:,2),...
        dis_fit_3guess(:,3),dis_fit_3guess(:,4),...
        linebroadening, zeropadsize);
    if(constrainRBCpeak)
        dis3peakFit = dis3peakFit.setBounds([0.1*rbcArea 0.1*barrierArea 0.1*gasArea],[inf inf inf],...
            [dis5peakFit.freq(1)-0.001 -inf -inf],[dis5peakFit.freq(1)+0.001 dis5peakFit.freq(1)-0.002 dis5peakFit.freq(1)-0.002],...
            [0 0 0],[400 400 400],...
            [-inf -inf -inf],[inf inf inf]);
    end
    dis3peakFit = dis3peakFit.fitTimeDomainSignal();
    
    % Store fit
    teFits_3peak{iTE} = dis3peakFit;
    
    % Check if fits are good
    goodFit(iTE) = checkFits(dis5peakFit, dis3peakFit, ref_freq);
end

%% Calculate TE90 if its not given from Dixon scan
te90s = zeros(nTE,1);
if(isempty(te_rbcBarrier))
    te_rbcBarrier = 0;
    totalGoodFits = 0;
    for iTE = 1:nTE
        if(goodFit(iTE))
            totalGoodFits = totalGoodFits + 1;
            
            % Calculate TE 90
            deltaF = teFits_3peak{iTE}.freq(2)-teFits_3peak{iTE}.freq(1);
            deltaPhase = teFits_3peak{iTE}.phase(2)-teFits_3peak{iTE}.phase(1);
            time180 = abs(1/(2*deltaF));
            te90s(iTE) = (90-deltaPhase)/(360*deltaF);
            while(te90s(iTE)>(minte+time180))
                % This te is too high, so subtract 180 deg of phase
                te90s(iTE) = te90s(iTE) - time180;
            end
            while(te90s(iTE)<minte)
                % This TE is too low, so add 180 deg of phase
                te90s(iTE) = te90s(iTE) + time180;
            end
            thisTE = te90s(iTE)*1E6
            te_rbcBarrier = te_rbcBarrier + te90s(iTE);
        end
    end
    te_rbcBarrier = te_rbcBarrier/totalGoodFits;
end

%% Convert fits to ppm units
for iTE=1:nTE
    % Convert 3-peak fits
    teFits_3peak{iTE}.freq = (teFits_3peak{iTE}.freq - ref_freq)/xeGam;
    teFits_3peak{iTE}.fwhm = teFits_3peak{iTE}.fwhm/xeGam;
    
    % Convert 5 peak fits
    teFits_5peak{iTE}.freq = (teFits_5peak{iTE}.freq - ref_freq)/xeGam;
    teFits_5peak{iTE}.fwhm = teFits_5peak{iTE}.fwhm/xeGam;
end

%% Display healthy reference spectra
specFitFig = figure('Name','Spectral Fitting','Units','Pixels','Position',[66 372 1200 600]);

% Calculate zeropadded time and frequencies
dwell_time = teFits_3peak{1}.t(2)-teFits_3peak{1}.t(1);
zeroPaddedTime = min(teFits_3peak{1}.t(:)) + dwell_time*((1:teFits_3peak{1}.zeroPadSize)-1)';
zeroPaddedFreq = linspace(-0.5,0.5,teFits_3peak{1}.zeroPadSize+1)/dwell_time;
zeroPaddedFreq = zeroPaddedFreq(1:(end-1)); % Take off last sample to have nSamples
zeroPaddedFreq = zeroPaddedFreq(:); % Make a column vector
zeroPaddedFreq_ppm = (zeroPaddedFreq - ref_freq)/xeGam; % Put into ppm units

% 3 peak fits
ax2 = subplot(2,1,2);
healthyRef_3peaks = NMR_Mix([0.271380049	0.480061072	0.059299147],...
    [216.182746	197.4973072	-0.208983034],...
    [10.00489504	8.609091632	2.111489336],...
    [0 0 0]);
individualSpectrums_ref3 = healthyRef_3peaks.calcComponentLorentzianCurves(zeroPaddedFreq_ppm)/(healthyRef_3peaks.area(2)/(pi*healthyRef_3peaks.fwhm(2)));
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_ref3(:,1)),'Linewidth',4,'color','k','LineStyle','-');
hold on;
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_ref3(:,2)),'Linewidth',4,'color','k','LineStyle','-');
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_ref3(:,3)),'Linewidth',4,'color','k','LineStyle','-');

% 5 peak fit
ax1 = subplot(2,1,1);
healthyRef_5peaks = NMR_Mix([0.311088817	0.237421553	0.451489631	0.042256641	0.036092003],...
    [216.6355279	201.5047859	196.5347491	0.175732399	-2.943695425],...
    [11.34666752	14.78528854	9.226784674	2.066920757	3.923064434],...
    [0 0 0 0 0]);
individualSpectrums_ref5 = healthyRef_5peaks.calcComponentLorentzianCurves(zeroPaddedFreq_ppm)/...
    ((healthyRef_5peaks.area(2)/(pi*healthyRef_5peaks.fwhm(2)))+...
    (healthyRef_5peaks.area(3)/(pi*healthyRef_5peaks.fwhm(3))));
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_ref5(:,1)),'Linewidth',4,'color','k','LineStyle','-');
hold on;
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_ref5(:,2)),'Linewidth',4,'color','k','LineStyle','-');
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_ref5(:,3)),'Linewidth',4,'color','k','LineStyle','-');
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_ref5(:,4)),'Linewidth',4,'color','k','LineStyle','-');
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_ref5(:,5)),'Linewidth',4,'color','k','LineStyle','-');

%% Calculate Spectrum from "good fits"
barrier_norm_ref3 = 0;
barrier_norm_ref5 = 0;
totalGoodFits = 0;
meanSpectrum_3peak = NMR_Mix([0 0 0],[0 0 0],[0 0 0],[0 0 0]);
meanSpectrum_5peak = NMR_Mix([0 0 0 0 0],[0 0 0 0 0],[0 0 0 0 0],[0 0 0 0 0]);
for iTE = 1:nTE
    if(goodFit(iTE))
        totalGoodFits = totalGoodFits + 1;
        
        meanSpectrum_3peak.area = meanSpectrum_3peak.area + teFits_3peak{iTE}.area;
        meanSpectrum_3peak.freq = meanSpectrum_3peak.freq + teFits_3peak{iTE}.freq;
        meanSpectrum_3peak.fwhm = meanSpectrum_3peak.fwhm + teFits_3peak{iTE}.fwhm;
        barrier_norm_ref3 = barrier_norm_ref3 + ...
            teFits_3peak{iTE}.area(2)/(pi*teFits_3peak{iTE}.fwhm(2));
        
        meanSpectrum_5peak.area = meanSpectrum_5peak.area + teFits_5peak{iTE}.area;
        meanSpectrum_5peak.freq = meanSpectrum_5peak.freq + teFits_5peak{iTE}.freq;
        meanSpectrum_5peak.fwhm = meanSpectrum_5peak.fwhm + teFits_5peak{iTE}.fwhm;
        barrier_norm_ref5 = barrier_norm_ref5 + ...
            teFits_5peak{iTE}.area(2)/(pi*teFits_5peak{iTE}.fwhm(2)) + ...
            teFits_5peak{iTE}.area(3)/(pi*teFits_5peak{iTE}.fwhm(3));
    end
end
meanSpectrum_3peak.area = meanSpectrum_3peak.area/barrier_norm_ref3;
meanSpectrum_3peak.freq = meanSpectrum_3peak.freq/totalGoodFits;
meanSpectrum_3peak.fwhm = meanSpectrum_3peak.fwhm/totalGoodFits;
meanSpectrum_5peak.area = meanSpectrum_5peak.area/barrier_norm_ref5;
meanSpectrum_5peak.freq = meanSpectrum_5peak.freq/totalGoodFits;
meanSpectrum_5peak.fwhm = meanSpectrum_5peak.fwhm/totalGoodFits;

%% Display average of "good" 3- and 5-peak fits
linestyles = {'-','--','-.',':'};

% Calculate 3-peak spectrum and display
individualSpectrums_3peak = meanSpectrum_3peak.calcComponentLorentzianCurves(zeroPaddedFreq_ppm);
ax2 = subplot(2,1,2);
hold on
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_3peak(:,1)),'Linewidth',3,'color',[0.8500    0.3250    0.0980],'LineStyle','-');
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_3peak(:,2)),'Linewidth',3,'color',[0.4660    0.6740    0.1880],'LineStyle','-');
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_3peak(:,3)),'Linewidth',3,'color',[0    0.4470    0.7410],'LineStyle','-');
ylabel('Component Intensity');
set(ax2,'XDir','reverse');

% Calculate 5-peak spectrum and display
individualSpectrums_5peak = meanSpectrum_5peak.calcComponentLorentzianCurves(zeroPaddedFreq_ppm);
ax1 = subplot(2,1,1);
plot(zeroPaddedFreq_ppm ,real(individualSpectrums_5peak(:,1)),'Linewidth',3,'color',[0.8500    0.3250    0.0980],'LineStyle','-');
hold on
plot(zeroPaddedFreq_ppm,real(individualSpectrums_5peak(:,2)),'Linewidth',3,'color',[ 0.9290    0.6940    0.1250],'LineStyle','-');
plot(zeroPaddedFreq_ppm,real(individualSpectrums_5peak(:,3)),'Linewidth',3,'color',[0.4660    0.6740    0.1880],'LineStyle','-');
plot(zeroPaddedFreq_ppm,real(individualSpectrums_5peak(:,4)),'Linewidth',3,'color',[0.3010    0.7450    0.9330],'LineStyle','-');
plot(zeroPaddedFreq_ppm,real(individualSpectrums_5peak(:,5)),'Linewidth',3,'color',[0    0.4470    0.7410],'LineStyle','-');
ylabel('Component Intensity');
set(ax1,'XDir','reverse');
linkaxes([ax1 ax2],'x');
xlim([-50 250]);
subplot(2,1,1);
ylim([ 0 1]);
subplot(2,1,2);
ylim([ 0 1.2]);

%% List out fit parameters
% 5-Peak
for iTE = 1:nTE
    % First normalize areas to barrier signal
    teFits_5peak{iTE}.area = teFits_5peak{iTE}.area/sum(teFits_5peak{iTE}.area(2:3));
    
    fprintf (['\n*** 5 Peak fit for TE = ' num2str(1E6*TEs(iTE)) ' usec ***\n']);
    displayPPMFits(teFits_5peak{iTE});
end
% 3-Peak
for iTE = 1:nTE
    % First normalize areas to barrier signal
    teFits_3peak{iTE}.area = teFits_3peak{iTE}.area/teFits_3peak{iTE}.area(2);
    
    fprintf (['\n*** 3 Peak fit for TE = ' num2str(1E6*TEs(iTE)) ' usec ***\n']);
    if(constrainRBCpeak)
        disp(['(RBC frequency constrained at ' num2str(teFits_5peak{iTE}.freq(1)) ' ppm)']);
    end
    displayPPMFits(teFits_3peak{iTE});
end

% Calculate RBC:barrier at various TEs
compSignal_3peak_875 = meanSpectrum_3peak.calcComponentTimeDomainSignal(875E-6);
compSignal_3peak_te90 = meanSpectrum_3peak.calcComponentTimeDomainSignal(te_rbcBarrier);
compSignal_5peak_875 = meanSpectrum_5peak.calcComponentTimeDomainSignal(875E-6);
compSignal_5peak_te90 = meanSpectrum_5peak.calcComponentTimeDomainSignal(te_rbcBarrier);
rbcBarrier_3peak_875 = abs(compSignal_3peak_875(1))/abs(compSignal_3peak_875(2));
rbcBarrier_3peak_te90 = abs(compSignal_3peak_te90(1))/abs(compSignal_3peak_te90(2));
rbcBarrier_5peak_875 = abs(compSignal_5peak_875(1))/(abs(compSignal_5peak_875(2))+abs(compSignal_5peak_875(2)));
rbcBarrier_5peak_te90 = abs(compSignal_5peak_te90(1))/(abs(compSignal_5peak_te90(2))+abs(compSignal_5peak_te90(3)));

%% Display Summary
fprintf ('\n*** SUMMARY ***\n');
if(any(~goodFit))
    warning(['Ignored data from echos :' num2str(find(~goodFit))]);
end
disp(['mean RBC freq = ' num2str(meanSpectrum_3peak.freq(1)) ' ppm']);
disp(['mean TE90 = ' num2str(te_rbcBarrier*1E6) ' usec']);
disp(['mean 3-peak RBC:Barrier(875usec) = ' num2str(mean(rbcBarrier_3peak_875))]);
disp(['mean 5-peak RBC:Barrier(875usec) = ' num2str(mean(rbcBarrier_5peak_875))]);
disp(['mean 3-peak RBC:Barrier(TE90) = ' num2str(mean(rbcBarrier_3peak_te90))]);
disp(['mean 5-peak RBC:Barrier(TE90) = ' num2str(mean(rbcBarrier_5peak_te90))]);
disp(['Flip angle ~' num2str(flip_angle) ' (' num2str(flip_err) ' error)']);


    function displayPPMFits(theFitObject)
        disp('area (arbs)  Freq (ppm)  Linewidth(ppm)  Phase(degrees)');
        for iComp = 1:length(theFitObject.area)
            disp([sprintf('%8.2f%%',100*theFitObject.area(iComp)) ' ' ...
                sprintf('%+8.2f',theFitObject.freq(iComp))  '  ' ...
                sprintf('%8.2f',theFitObject.fwhm(iComp)-theFitObject.lineBroadening) '  ' ...
                sprintf('%+9.2f',theFitObject.phase(iComp))]);
        end
    end

    function isGoodFit = checkFits(fits_5peak, fits_3peak, freq_ref)
        % Assume its good, unless proven otherwise
        isGoodFit = 1;
        
        if(any(fits_3peak.area < 0) || any(fits_5peak.area < 0))
            isGoodFit = 0;
            return;
        end
        
        % Frequencies
        if(((fits_3peak.freq(1)-freq_ref)/xeGam > 219) || ...
                ((fits_3peak.freq(1)-freq_ref)/xeGam < 213) || ...
                ((fits_5peak.freq(1)-freq_ref)/xeGam > 219) || ...
                ((fits_5peak.freq(1)-freq_ref)/xeGam < 213))
        	isGoodFit=0;
            return
        end
        if(((fits_3peak.freq(2)-freq_ref)/xeGam > 205) || ...
                ((fits_3peak.freq(2)-freq_ref)/xeGam < 180) || ...
                ((fits_5peak.freq(2)-freq_ref)/xeGam > 215) || ...
                ((fits_5peak.freq(2)-freq_ref)/xeGam < 198) || ...
                ((fits_5peak.freq(3)-freq_ref)/xeGam > 198) || ...
                ((fits_5peak.freq(3)-freq_ref)/xeGam < 180))
            isGoodFit=0;
            return
        end
        if(((fits_3peak.freq(3)-freq_ref)/xeGam > 15) || ...
                ((fits_3peak.freq(3)-freq_ref)/xeGam < -15) || ...
                ((fits_5peak.freq(4)-freq_ref)/xeGam > 15) || ...
                ((fits_5peak.freq(4)-freq_ref)/xeGam < -15) || ...
                ((fits_5peak.freq(5)-freq_ref)/xeGam > 15) || ...
                ((fits_5peak.freq(5)-freq_ref)/xeGam < -15))
            isGoodFit=0;
            return
        end
        
        % FWHMs
        if(any(fits_3peak.fwhm < 0) || any(fits_5peak.fwhm < 0) ||...
                any(fits_3peak.fwhm/xeGam > 25) || any(fits_5peak.fwhm/xeGam > 25))
            isGoodFit = 0;
            return;
        end
    end
end