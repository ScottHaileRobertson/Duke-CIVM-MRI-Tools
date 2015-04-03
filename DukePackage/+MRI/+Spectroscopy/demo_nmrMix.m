% clc; close all; clear all;
dwell_time = 62E-6; %sec
npts = 512;
t = (dwell_time*((1:npts) - 1)); %sec

tau = [5E-3 3E-3 12E-3]';    %(sec)T2*=2msec
amp = [4 2 5]';    %arbs
freq = [-4E3 30 400 ]';   %Hz
phase = [0 -20 45]';     %deg
fwhm = 1./(pi*tau) %Hz
nmrMix = NMR_Mix(amp, freq, fwhm, phase);

% Calculate temporal frequencies (for FFT)
f = NMR_Mix.calcFftFreq(t);

% Calculate signal and spectrum
signal = sum(nmrMix.calcTimeDomainSignal(t),2);
spectrum = sum(nmrMix.calcSpectralDomainSignal(f),2);

% add noise
noisySignal = signal + 7*([rand(size(signal)) + 1i*rand(size(signal))]-[0.5+0.5i]);

% Fit data in time domain
fitMix = NMR_Mix([],[],[],[]);
fitMix = fitMix.autoFitComponents(noisySignal, t, 3);
timefit_signal = sum(fitMix.calcTimeDomainSignal(t),2);
timefit_spectrum = sum(fitMix.calcSpectralDomainSignal(f),2);

% % Fit data in spectral domain
% [spectralfit_amp spectralfit_freq spectralfit_fwhm spectralfit_phase] = ...
% 	fitSpectralDomainSignal(spectrum, f, guess_amp, guess_freq, guess_fwhm, guess_phase)
% spectralfit_signal = calcTimeDomainSignal(spectralfit_amp,spectralfit_freq,spectralfit_fwhm,spectralfit_phase,t);
% spectralfit_signal = sum(spectralfit_signal,2);
% spectralfit_spectrum = calcSpectralDomainSignal(spectralfit_amp,spectralfit_freq,spectralfit_fwhm,spectralfit_phase,f);
% spectralfit_spectrum = sum(spectralfit_spectrum,2);

figure();
subplot(4,1,1);
plot(t, real(noisySignal),'-b');
hold on
plot(t,real(timefit_signal),'--g');
% plot(t,real(spectralfit_signal),':r');
xlabel('Time (sec)');
ylabel('Intensity')
title('Real Time Domain Signal');
legend('Ideal','Time fit','Spectral fit');

subplot(4,1,2);
plot(t, imag(noisySignal),'-b');
hold on
plot(t,imag(timefit_signal),'--g');
% plot(t,imag(spectralfit_signal),':r');
xlabel('Time (sec)');
ylabel('Intensity')
title('Imaginary Time Domain Signal');
legend('Ideal','Time fit','Spectral fit');

subplot(4,1,3);
plot(f, real(spectrum),'-b');
hold on
plot(f,real(timefit_spectrum),'--g');
% plot(f,real(spectralfit_spectrum),':r');
xlabel('Spectral Frequency (Hz)');
ylabel('Intensity')
title('Real Spectral Domain Signal');
legend('Ideal','Time fit','Spectral fit');

subplot(4,1,4);
plot(f, imag(spectrum),'-b');
hold on
plot(f,imag(timefit_spectrum),'--g');
% plot(f,imag(spectralfit_spectrum),':r');
xlabel('Spectral Frequency (Hz)');
ylabel('Intensity')
title('Imaginary Time Domain Signal');
legend('Ideal','Time fit','Spectral fit');