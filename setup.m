% Get current path
rootDistribDir = pwd();

% Add Recon package to path
disp('Adding Duke package to MATLAB path...');
path(genpath([rootDistribDir filesep() 'DukePackage']),path);
