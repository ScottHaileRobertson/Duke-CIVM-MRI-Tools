% Find pfile
pfile_path = filepath('C:\Users\Scott\Desktop\rohan_20150331\P33792.7')

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

% Create array of sample times (sec)
npts = pfile.rdb.rdb_hdr_frame_size;                   % Number of samples
bw = 1000*pfile.rdb.rdb_hdr_user12;                         % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw);                                 % Time between each sample
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
t = dwell_time*(0:(npts-1));

% Fit spectra
nmrMix = NMR_Mix([1921 496 300],[-59.5 -20.7 -4600],[131.44 30 40],[-45 25 0],pfile.rdb.rdb_hdr_ps_mps_freq/10);
nmrMix = nmrMix.fitTool(pfile.data, t)
