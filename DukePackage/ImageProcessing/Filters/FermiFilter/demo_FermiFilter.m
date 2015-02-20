clc; close all;

load('mri'); % Matlab default
D = double(squeeze(D));

% Show volume before filtration
figure();
imslice(abs(D),'Unfiltered volume');

% Fermi filter parameters
fermi_rolloff=0.15;
fermi_radius=0.75;

% Apply fermi filter
b = FermiFilter(D, fermi_rolloff, fermi_radius);

% Show filtered volume
figure();
imslice(abs(b),'Filtered volume');
