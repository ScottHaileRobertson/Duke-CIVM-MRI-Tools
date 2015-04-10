
% pfile_path = filepath('/home/scott/Public/data/20150311')
% pfile_path = '/home/scott/Desktop/Dixon/P12800.7';
% pfile_path = filepath('C:\Users\Scott\Desktop\subject002_065\P32768.7');
pfile_path = filepath('C:\Users\Scott\Desktop\P06144.7');
% pfile_path = filepath('C:\Users\Scott\');

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

    zeropadsize = 2048*4;
    linebroadening = 0; %Hz

% Split disolved data into separate pfile
nDisFrames = 200;
nGasFrames = 1;
dis_pfile = pfile;
dis_pfile.data = pfile.data(:,1:nDisFrames);
dis_pfile.rdb.rdb_hdr_user20 = nDisFrames; % nframes
center_freq = pfile.rdb.rdb_hdr_ps_mps_freq/10;

% Split dissolved data into 4 pfiles (one for each TE)
nTE = 4;
TEs = [875 975 1075 1175];
disTE_pfile = cell(1,nTE);
b = zeros(4,2048);
for iTE = 1:nTE
	disTE_pfile{iTE} = dis_pfile;
	disTE_pfile{iTE}.data = disTE_pfile{iTE}.data(:,iTE:nTE:end);
	disTE_pfile{iTE}.rdb.rdb_hdr_user20 = size(disTE_pfile{iTE}.data,2); % nframes
	npts = disTE_pfile{iTE}.rdb.rdb_hdr_frame_size;
	
	% Average
% 	avgdata = disTE_pfile{iTE}.data(:,1);
	avgdata = mean(disTE_pfile{iTE}.data,2);
	
% 	% calculate time domain baselines
% 	base_points=round(baseline_fraction*npts);
% 	i_baseline=mean(imag(avgdata(end-base_points:end)));
% 	r_baseline=mean(real(avgdata(end-base_points:end)));
% 	avgdata = avgdata-(r_baseline+i*i_baseline);      % subtract baselines
	
	% Pull relavent info from header
	bw = pfile.rdb.rdb_hdr_user12;
	npts = pfile.rdb.rdb_hdr_frame_size;
	
% 	% Line broaden
	dwell_time = 62E-6;%1/(2*bw*1000);
    broaddata = avgdata;
% 	t = (TEs(iTE)*1E-6)+(dwell_time*((1:npts) - 1)); %sec
    t = (dwell_time*((1:npts) - 1)); %sec
    
    nComp = 3;
    
%     %            Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
%     fit_guess = [   1           -7.77           175          0; % Component #1
% %                     1           -7           100          0; % Component #2
% %                     1            -3832           100          0; % Component #3
%                     1            -348           100          0; % Component #4
%                     1           -3843           100          0]; % Component #5
% %                     1          -3992           100          0]; % Component #
%                 
%     nmrMix = NMR_Mix(broaddata, t, fit_guess(:,1),fit_guess(:,2),...
%         fit_guess(:,3),fit_guess(:,4),zeropadsize,linebroadening,center_freq);
    nmrMix = NMR_Mix(broaddata, t,[],[],[],[],zeropadsize, linebroadening, center_freq);
    
%     nmrMix = nmrMix.fitTimeDomainSignal();
%     nmrMix = nmrMix.fitTool();
    nmrMix = nmrMix.autoAddComponents(nComp);
    
    figure();
    nmrMix.plotFit();
    title(['TE' num2str(iTE)]);
    
%     % Using guess from averaged data, fit each spectrum
%     nSpec = size(disTE_pfile{iTE}.data,2);
%     freqs = zeros(nComp,nSpec);
%     fwhms = zeros(nComp,nSpec);
%     phases = zeros(nComp,nSpec);
%     for iSpec = 1:nSpec
%         indivData = disTE_pfile{iTE}.data(:,iSpec);
%         thisMix = nmrMix;
%         thisMix = thisMix.fitTimeDomainSignal(indivData,t);
%         freqs(:,iSpec) = thisMix.freq;
%         fwhms(:,iSpec) = thisMix.fwhm;
%         phases(:,iSpec) = thisMix.phase;
%     end
%     
%     figure(20);
%     subplot(5,4,0+iTE);
%     plot(repmat(1:nSpec,[nComp 1])',freqs');
%     xlabel('Frame number');
%     ylabel('Frequency (Hz)')
%     title(['TE' num2str(iTE)]);
%     
%     subplot(5,4,4+iTE);
%     plot(repmat(1:nSpec,[nComp 1])',(freqs - repmat(freqs(:,1),[1 nSpec]))');
%     xlabel('Frame number');
%     ylabel('Differential Frequency (Hz)')
%     
%     subplot(5,4,8+iTE);
%     plot(repmat(1:nSpec,[nComp 1])',fwhms');
%     xlabel('Frame number');
%     ylabel('FWHM (Hz)')
%     
%         subplot(5,4,12+iTE);
%     plot(repmat(1:nSpec,[nComp 1])',(fwhms - repmat(fwhms(:,1),[1 nSpec]))');
%     xlabel('Frame number');
%     ylabel('Differential FWHM (Hz)')
%     
%     subplot(5,4,16+iTE);
%     plot(repmat(1:nSpec,[nComp 1])',phases');
%     xlabel('Frame number');
%     ylabel('Phases (deg)')
%     nmrMix
    nmrMixes{iTE} = nmrMix;
end

% Calculate ratios
colors = get(groot,'defaultAxesColorOrder');
area_rbc = nmrMixes{iTE}.amp.*nmrMixes{iTE}.fwhm;
area_barrier = nmrMixes{iTE}.amp.*nmrMixes{iTE}.fwhm;
area_gas = nmrMixes{iTE}.amp.*nmrMixes{iTE}.fwhm;
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
plot(repmat(TEs,[3 1])',[area_rbc area_barrier area_gas]);
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
gasMix = gasMix.fitTimeDomainSignal()
% gasMix = NMR_Mix(gas_pfile.data,t,[],[],[],[],zeropadsize,linebroadening,center_freq);
% gasMix = gasMix.fitTool();
figure()
gasMix.plotFit();

% Split flip calibration frames into separate pfile
flipCal_pfile = pfile;
flipCal_pfile.data = pfile.data(:,nDisFrames+nGasFrames+1:end);
nFlipCal = size(flipCal_pfile.data,2);
% flipCal_pfile.rdb_hdr_user20 = nFlipCal; % nframes





% pfile.data = linebroaden(pfile.data, linebroadening = 50)



% Pick a reasonable DC sample
% dc_sample_idx = MRI.DataProcessing.calcDCsample(max(radialDistance,[],2));

% % Calculate flip angle
% MRI.DataProcessing.calcFlipAngle(pfile, dc_sample_idx);

