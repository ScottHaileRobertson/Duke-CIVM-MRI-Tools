% This is a demo script that will fit a spectroscopy file 

% Find pfile
pfile_path = filepath('C:\Users\Scott\Desktop\rohan_20150331\')

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(pfile_path);
displayPfileHeaderInfo(pfile);

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
pfile = MRI.DataProcessing.removeBaselineViews(pfile);

% discard first view, update header
pfile.data = pfile.data(:,2:end);
pfile.rdb.rdb_hdr_user20 = size(pfile.data,2);

% Average all frames
pfile.data = mean(pfile.data,2);
pfile.rdb.rdb_hdr_user20 = 1; % Update header to only have one averaged frame

%% Prepare for spectroscopy fitting
%Create array of sample times (sec)
npts = pfile.rdb.rdb_hdr_frame_size;                   % Number of samples
bw = 1000*pfile.rdb.rdb_hdr_user12;                    % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw);                                 % Time between each sample
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
t = dwell_time*(0:(npts-1));

%% Fit spectra
%            Amplitude   Frequency(Hz)   FWHM(Hz)    Phase(deg)
fit_guess = [   1             102           100          0; % Component #1
                1            31           100          0; % Component #2
                1           -8           100          0; % Component #3
                1           -75           100          0;% Component #4
                1           -110           100          0;% Component #5
                1           -4646           100          0;% Component #6
                1           -4725           100          0];% Component #7
center_freq = pfile.rdb.rdb_hdr_ps_mps_freq/10;
nmrMix = NMR_Mix(fit_guess(:,1),fit_guess(:,2),fit_guess(:,3),fit_guess(:,4),center_freq);
% nmrMix = NMR_Mix([],[],[],[],center_freq);

% Option one: just fit the guess peaks
nmrMix = nmrMix.fitTimeDomainSignal(pfile.data,t);

% % Option two: use the fit tool
% nmrMix = nmrMix.fitTool(pfile.data, t)

%% Display the final fit
nmrMix.displayFit(pfile.data,t);

