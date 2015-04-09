% clc; close all; clear all;
dwell_time = 62E-6; %sec
npts = 512;
t = (dwell_time*((1:npts) - 1)); %sec

% tau = 1./(pi*[5E-3 3E-3 1E-3])';    %(sec)T2*=2msec
amp = [10 10 10]';    %arbs
freq = [-4E3 30 400 ]';   %Hz
phase = [0 -20 45]';     %deg
% fwhm = 1./(pi*tau) %Hz
fwhm = [100 200 300];
nmrMix = NMR_Mix(amp, freq, fwhm, phase);

% Calculate temporal frequencies (for FFT)
f = NMR_Mix.calcFftFreq(t);

% Calculate signal and spectrum
signal = sum(nmrMix.calcTimeDomainSignal(t),2);
spectrum = sum(nmrMix.calcSpectralDomainSignal(f),2);

% add noise
noisySignal = signal;% + 0*([rand(size(signal)) + 1i*rand(size(signal))]-[0.5+0.5i]);

% Fit data in time domain
fitMix = NMR_Mix([],[],[],[]);
fitMix = fitMix.autoAddComponents(noisySignal, t, 3);

freq_fit = fitMix.freq
phase_fit = fitMix.phase
amp = fitMix.amp
fwhm = fitMix.fwhm

% % Fit data in spectral domain
% [spectralfit_amp spectralfit_freq spectralfit_fwhm spectralfit_phase] = ...
% 	fitSpectralDomainSignal(spectrum, f, guess_amp, guess_freq, guess_fwhm, guess_phase)
% spectralfit_signal = calcTimeDomainSignal(spectralfit_amp,spectralfit_freq,spectralfit_fwhm,spectralfit_phase,t);
% spectralfit_signal = sum(spectralfit_signal,2);
% spectralfit_spectrum = calcSpectralDomainSignal(spectralfit_amp,spectralfit_freq,spectralfit_fwhm,spectralfit_phase,f);
% spectralfit_spectrum = sum(spectralfit_spectrum,2);

fitMix.plotFit(noisySignal,t);