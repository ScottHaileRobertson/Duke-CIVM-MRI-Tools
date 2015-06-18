%% NMR_MIX
% A class that represents a series of exponentially decaying components. 
% The class knows how to calculate time or spectral domain signals. The 
% key assumption is that all components experience perfectly exponential 
% decay. This can be modeled as
%
% s(t) = amp*exp(i*phase)*exp(-pi*fwhm*t)*exp(-i*2*pi*freq*t)
% S(f) = (amp*exp(i*phase))/(i*2*pi*(freq-f)-pi*fwhm)
%
%
classdef NMR_Mix < handle
    
    properties
        amp;   % Amplitude of each component (arbs)
        freq;  % Frequency of each component in Hz
        fwhm;  % FWHM of each component in Hz
        phase; % Phase of each component in degrees
        
    end
    
    methods
        function obj = NMR_Mix(amp, freq, fwhm, phase)
            % Constructs an NMR_Mix object with the given time domain
            % signal  and optional knowledge of component amplitudes
            % (arbs), frequencies (in Hz), fwhms (in Hz), and phases (in
            % radians). 
            obj = obj.addComponents(amp, freq, fwhm, phase);
        end
        
        function obj = resetComponents(obj,amplitude, frequency, fwhm, phase)
            % Resets component(s) to NMR_Mix object, then resorts to keep
            % frequencies in order
            obj.amp = amplitude(:)';
            obj.freq = frequency(:)';
            obj.fwhm = fwhm(:)';
            obj.phase = phase(:)';
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = addComponents(obj,amplitude, frequency, fwhm, phase)
            % Adds component(s) to NMR_Mix object, then resorts to keep
            % frequencies in order
            obj.amp = [obj.amp amplitude(:)'];
            obj.freq = [obj.freq frequency(:)'];
            obj.fwhm = [obj.fwhm fwhm(:)'];
            obj.phase = [obj.phase phase(:)'];
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = removeComponents(obj, componentNumbers)
            % Removes the component(s) given their index number
            obj.amp(componentNumbers) = [];
            obj.freq(componentNumbers) = [];
            obj.fwhm(componentNumbers) = [];
            obj.phase(componentNumbers) = [];
        end
        
        function obj = sortByFreq(obj)
            % Sorts all components by thier frequencies in decending order
            [sortFreq sortIdx] = sort(obj.freq, 'descend');
            obj.amp = obj.amp(sortIdx);
            obj.freq = obj.freq(sortIdx);
            obj.fwhm = obj.fwhm(sortIdx);
            obj.phase = obj.phase(sortIdx);
        end
        
        function componentTimeDomainSignal = calcComponentTimeDomainSignal(obj,t)
            % Calculates the time domain signal from the individual components
            % of the NMR mix at the given time points (t is in sec). Note, 
            % this function returns the time signal for each individual
            % component. The overall "mix" signal can be obtained with 
            % calcTimeDomainSignal
            t = t(:); % immune to passing in a row or col vec
            nt = length(t);
            componentTimeDomainSignal = repmat(obj.amp,[nt 1]).*...
                exp(1i*pi/180*repmat(obj.phase,[nt 1]) - ...
                pi*t*obj.fwhm + 1i*2*pi*t*obj.freq);
        end
        
        function timeDomainSignal = calcTimeDomainSignal(obj,t)
            % Calculates the net time domain signal from the mix of NMR
            % components at the given time points (t is in sec). Note, this
            % function returns the net time signal from the sum of all 
            % components in the NMR_Mix. The individual component signals
            % can be obtained with calcComponentTimeDomainSignal
            timeDomainSignal = sum(obj.calcComponentTimeDomainSignal(t),2);
        end
        
        function componentSpectralDomainSignal = calcComponentSpectralDomainSignal(obj,f)
            % Calculates the spectral domain signal from the individual 
            % components of the NMR_mix at the given spectral frequencies 
            % (f is in Hz). Note, this function returns the spectrum for each individual
            % component. The overall "mix" spectrum can be obtained with 
            % calcSpectralDomainSignal.
            f = f(:); % immune to passing in a row or col vec
            nf = length(f);
            componentSpectralDomainSignal = (repmat(obj.amp,[nf 1]).*exp(1i*pi/180*repmat(obj.phase,[nf 1])))./...
                (1i*2*pi*(repmat(f,[1 length(obj.amp)])-repmat(obj.freq,[nf 1]))+...
                pi*repmat(obj.fwhm,[nf 1]));
        end
        
        function spectralDomainSignal = calcSpectralDomainSignal(obj,f)
            % Calculates the spectral domain signal from the mix of NMR 
            % components at the given spectral frequencies (f is in Hz). 
            % Note, this function returns the net spectrum from the sum of
            % all components in the NMR_Mix. The individual component 
            % spectrums can be obtained with calcComponentSpectralDomainSignal.
            spectralDomainSignal = sum(obj.calcComponentSpectralDomainSignal(f),2);
        end
        
        function componentLorentzianCurves = calcComponentLorentzianCurves(obj,f)
            % Calculates unphased lorenzian curves for each NMR component 
            % at the given spectral frequencies (f is in Hz). The 
            % lorentzians are just the real part of the spectrum.
            % Note, this function returns a curve for each individual
            % component. The overall "mix" of lorentzian curves can be
            % obtained with calcLorentzianCurves.
            f = f(:); % immune to passing in a row or col vec
            nf = length(f);
            nc = length(obj.amp);
            linewidth = pi*repmat(obj.fwhm,[nf 1]);
            componentLorentzianCurves = (repmat(obj.amp,[nf 1]).*linewidth)./...
                ((2*pi*(repmat(obj.freq,[nf 1])-repmat(f,[1 nc]))).^2+linewidth.^2);
        end
        
        function lorentzianCurves = calcLorentzianCurves(obj,f)
            % Calculates unphased lorenzian curves for the net NMR_mix
            % at the given spectral frequencies (f is in Hz). The 
            % lorentzians are just the real part of the spectrum.
            % Note, this function returns a curve for the net signal from 
            % all components. The individual lorentzian curves can be
            % obtained with calcComponentLorentzianCurves.
            componentLorentzianCurves = sum(obj.calcComponentLorentzianCurves(f),2);
        end
        
        function componentDispersiveCurves = calcComponentDispersiveCurves(obj,f)
            % Calculates unphased dispersive curves for each NMR component 
            % at the given spectral frequencies (f is in Hz). The 
            % dispersive functions are just the imaginary part of the spectrum.
            % Note, this function returns a curve for each individual
            % component. The overall "mix" of dispersive functions can be
            % obtained with calcDispersiveCurves.
            f = f(:); % immune to passing in a row or col vec
            nf = length(f);
            nc = length(obj.amp);
            diff_freq = 2*pi*(repmat(obj.freq,[nf 1])-repmat(f,[1 nc]));
            componentDispersiveCurves = (repmat(obj.amp,[nf 1]).*diff_freq)./...
                (diff_freq.^2+(pi*repmat(obj.fwhm,[nf 1])).^2);
        end
        function dispersiveCurves = calcDispersiveCurves(obj,f)
            % Calculates unphased dispersive functions for the net NMR_mix
            % at the given spectral frequencies (f is in Hz). The 
            % dispersive functions are just the imaginary part of the spectrum.
            % Note, this function returns a curve for the net signal from 
            % all components. The individual dispersive functions can be
            % obtained with calcComponentDispersiveCurves.
            f = f(:); % immune to passing in a row or col vec
            nf = length(f);
            nc = length(obj.amp);
            diff_freq = 2*pi*(repmat(obj.freq,[nf 1])-repmat(f,[1 nc]));
            dispersiveCurves = (repmat(obj.amp,[nf 1]).*diff_freq)./...
                (diff_freq.^2+(pi*repmat(obj.fwhm,[nf 1])).^2);
        end
    end
end