
% pfile_path = filepath('/home/scott/Public/data/20150311')
% pfile_path = '/home/scott/Desktop/Dixon/P12800.7';
pfile_path = 'P12800.7';

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% Split disolved data into separate pfile
nDisFrames = 200;
nGasFrames = 1;
dis_pfile = pfile;
dis_pfile.data = pfile.data(:,1:nDisFrames);
dis_pfile.rdb.rdb_hdr_user20 = nDisFrames; % nframes

% Split dissolved data into 4 pfiles (one for each TE)
nTE = 4;
disTE_pfile = cell(1,nTE);
b = zeros(4,2048);
for iTE = 1:nTE
	disTE_pfile{iTE} = dis_pfile;
	disTE_pfile{iTE}.data = disTE_pfile{iTE}.data(:,iTE:nTE:end);
	disTE_pfile{iTE}.rdb.rdb_hdr_user20 = size(disTE_pfile{iTE}.data,2); % nframes
	npts = disTE_pfile{iTE}.rdb.rdb_hdr_frame_size;
	
	% Average
	avgdata = disTE_pfile{iTE}.data(:,1);
% 	avgdata = mean(disTE_pfile{iTE}.data,2);
	
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
% 	broaddata = linebroaden(disTE_pfile{iTE},linebroadening, avgdata);
    broaddata = avgdata;
    
	% Zeropad then fft
	fftdata = fftshift(fft(broaddata));
	t = (dwell_time*((1:npts) - 1)); %sec
    f = NMR_Mix.calcFftFreq(t);
    
    nmrMix = NMR_Mix([1 1 1],[30 300 4600],[50 50 50],[0 0 0],pfile.rdb.rdb_hdr_ps_mps_freq/10);
    nmrMix = nmrMix.fitTool(broaddata, t);
%     nmrMix = nmrMix.fitNewTimeDomainSignal(broaddata, t, 3);
    timefit_signal = sum(nmrMix.calcTimeDomainSignal(t),2);
    timefit_spectrum = sum(nmrMix.calcSpectralDomainSignal(f),2);
end
figure();
plot(b')
legend('TE=875us','TE=975us','TE=1075us','TE=1175us')

% Split Gas spectra into separate pfile
gas_pfile = pfile;
gas_pfile.data = pfile.data(:,nDisFrames + 1:nGasFrames);
gas_pfile.rdb.rdb_hdr_user20 = nGasFrames; % nframes

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

