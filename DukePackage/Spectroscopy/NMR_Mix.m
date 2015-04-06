%% NMR_MIX
% A class that fits or simulates a series of exponentially decating
% components. The class knows how to fit to a time domain signal, can
% calculate time or spectral domain signals, and can convert to various
% units. The key assumption is that all components are perfectly
% exponential. Thus ____
%
%

classdef NMR_Mix < matlab.mixin.CustomDisplay
    
    properties
        amp;   % Amplitude of each component (arbs)
        freq;  % Frequency of each component in Hz/Rad
        fwhm;  % FWHM of each component in Hz
        phase; % Phase of each component in radians
        fref = [];
    end
    
    methods
        function obj = NMR_Mix(amplitudes, frequencies, fwhms, phases, varargin)
            % Constructs an NMR_Mix object with the given amplitudes
            % (arbs),frequencies (in Hz), fwhms (in Hz), and phases (in
            % degrees). To start with an empty object (ex for fitting),
            % just use NMR_mix([],[],[],[])
            obj.amp = amplitudes(:);
            obj.freq = frequencies(:); % In Hz
            obj.fwhm = fwhms(:);       % In Hz
            obj.phase = phases(:);     % In Degrees
            if(nargin == 5)
                obj.fref = varargin{1};
            end
        end
        
        function obj = addComponent(obj,amplitude, frequency, fwhm, phase)
            obj.amp = [obj.amp; amplitude];
            obj.freq = [obj.freq; frequency];
            obj.fwhm = [obj.fwhm; fwhm];
            obj.phase = [obj.phase; phase];
        end
        
        function obj = removeComponent(obj, componentNumber)
            obj.amp(componentNumber) = [];
            obj.freq(componentNumber) = [];
            obj.fwhm(componentNumber) = [];
            obj.phase(componentNumber) = [];
        end
        
        function timeDomainSignal = calcTimeDomainSignal(obj,t)
            % Calculates the time domain signal from the mix of NMR
            % components at the given time points (t is in sec). Note, this
            % function returns the time signal for each individual
            % component. The overall "mix" signal can be obtained by
            % summing accross all components.
            
            % Convert to the matrix form
            compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            % Calculate time domain signal using matrix form
            timeDomainSignal= obj.calcTimeDomainSignalFromMatrix(compMat,t);
        end
        
        function spectralDomainSignal = calcSpectralDomainSignal(obj,f)
            % Calculates the spectral domain signal from the mix of NMR
            % components at the given spectral frequencies (f is in Hz).
            % Note, this function returns the spectrum for each individual
            % component. The overall "mix" spectrum can be obtained by
            % summing accross all components.
            
            % Convert to the matrix form
            compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            % Calculate spectral domain signal using matrix form
            spectralDomainSignal = obj.calcSpectralDomainSignalFromMatrix(compMat,f);
        end
        
        function lorentzianCurves = calcLorentzianCurves(obj,f)
            % Calculates unphased lorenzian curves from the mix of NMR
            % components at the given spectral frequencies (f is in Hz).
            % Note, this function returns a curve for each individual
            % component. The overall "mix" of lorentzian curves can be
            % obtained by summing accross all components.
            
            % convert to matrix form
            compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            lorentzianCurves = obj.calcLorentzianCurvesFromMatrix(compMat, f);
        end
        
        function dispersiveCurves = calcDispersiveCurves(obj,f)
            % Calculates unphased dispersive curves from the mix of NMR
            % components at the given spectral frequencies (f is in Hz).
            % Note, this function returns a curve for each individual
            % component. The overall "mix" of dispersive curves can be
            % obtained by summing accross all components.
            
            % convert to matrix form
            compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            dispersiveCurves = obj.calcDispersiveCurvesFromMatrix(compMat, f);
        end
        
        function [amp freq fwhm phase] = fitComponentToResidual(obj,timeDomainSignal,t)
            % Fits a single exponentially decaying component to the
            % residual signal. This is useful in peak fitting. Note,
            % however, that this function only fits based on residuals,
            % thus you are fitting the residuals, and not the signal. After
            % fitting the residuals, it is recommended that you "refit" all
            % the components as they will likely change slightly to better
            % accomodate the new compnent.
            
            % Calculate spectrum and spectral frequencies
            dwell_time = t(2)-t(1);
            measuredSpectrum = dwell_time*fftshift(fft(timeDomainSignal));
            f = obj.calcFftFreq(t);
            
            % If there are already components, fit them first
            if(~isempty(obj.amp) | ~isempty(obj.freq) | ...
                    ~isempty(obj.fwhm) | ~isempty(obj.phase))
                obj = obj.fitTimeDomainSignal(timeDomainSignal, t);
                currentSpectrum = measuredSpectrum - sum(obj.calcSpectralDomainSignal(f),2);
            else
                currentSpectrum = measuredSpectrum;
            end
            
            % Create a list of phased spectra
            nPhases = 360;
            nSamples = length(timeDomainSignal);
            phases = linspace(-pi,pi,nPhases+1);
            phases = phases(1:(end-1));
            phasedSpectraMatrix = repmat(currentSpectrum,[1 nPhases]);
            phasesMatrix = repmat(phases,[nSamples 1]);
            phasedSpectraMatrix = phasedSpectraMatrix.*exp(1i*phasesMatrix);
            
            % Find spectral frequency and rough phase for the largest peak
            % accross all phases - note (abs(real())) does a better job of
            % separating close peaks that are out of phase
            [maxVal maxIdx] = max(abs(phasedSpectraMatrix(:)));
            [maxFreqIdx maxPhaseIdx] = ind2sub([nSamples nPhases],maxIdx);
            maxFreq = f(maxFreqIdx);
            maxPhase = phases(maxPhaseIdx);
            
            % Use curve fitting to fit residual signal
            residualTimeDomainSignal = timeDomainSignal - sum(obj.calcTimeDomainSignal(t),2);
            largestPeakComponent = NMR_Mix(1,maxFreq,10,maxPhase); %Guess amplitude and FWHM for now
            largestPeakComponent = largestPeakComponent.fitTimeDomainSignal(residualTimeDomainSignal, t);
            
            % Return fit parameters
            amp   = largestPeakComponent.amp;
            freq  = largestPeakComponent.freq;
            fwhm  = largestPeakComponent.fwhm;
            phase = largestPeakComponent.phase;
        end
        
        function obj = autoAddComponent(obj,timeDomainSignal,t)
            % This function attemps to automatically fit and add a new
            % exponentially decaying signal to the NMR_mix object
            
            % Fit next component to residual signal, then add it to this
            % NMR_Mix object
            [add_amp add_freq, add_fwhm add_phase] = ...
                obj.fitComponentToResidual(timeDomainSignal,t);
            
            % Add fitted component to NMR_mix
            obj.amp = [obj.amp; add_amp];
            obj.freq = [obj.freq; add_freq];
            obj.fwhm = [obj.fwhm; add_fwhm];
            obj.phase = [obj.phase; add_phase];
            
            % Refit all components after addition of latest component
            obj = obj.fitTimeDomainSignal(timeDomainSignal, t);
        end
        
        function obj = autoAddComponents(obj, timeDomainSignal, t, nComponents)
            % This function attemps to automatically fit and add
            % n (defined by nComponents) new exponentially decaying signals
            % to the NMR_mix object
            
            nIter = nComponents - length(obj.amp);
            for iComp = 1:nIter
                obj = obj.autoAddComponent(timeDomainSignal,t);
            end
        end
        
        function obj = fitTool(obj,timeDomainSignal,t)
            % This function attemps to adds exponentially decaying signal
            % components until the user stops the fit.
            
            % Calculate spectrum
            dwell_time = t(2)-t(1);
            f = NMR_Mix.calcFftFreq(t);
            measuredSpectrum = dwell_time*fftshift(fft(timeDomainSignal));
            
            % Create a figure to show results
            thisFigure = figure();
            
            % Keep fitting until user says stop
            continueFitting = true;
            while(continueFitting)
                % Fit next component to residual signal, then add it to this
                % NMR_Mix object
                [add_amp add_freq, add_fwhm add_phase] = ...
                    obj.fitComponentToResidual(timeDomainSignal,t);
                
                % Create temporary NRM_mix object and add newest component
                temp_NMR_Mix = NMR_Mix(obj.amp, obj.freq, obj.fwhm, obj.phase);
                temp_NMR_Mix.amp = [temp_NMR_Mix.amp; add_amp];
                temp_NMR_Mix.freq = [temp_NMR_Mix.freq; add_freq];
                temp_NMR_Mix.fwhm = [temp_NMR_Mix.fwhm; add_fwhm];
                temp_NMR_Mix.phase = [temp_NMR_Mix.phase; add_phase];
                
                % Refit all components after addition of latest component
                temp_NMR_Mix = temp_NMR_Mix.fitTimeDomainSignal(timeDomainSignal, t);
                
                
                temp_NMR_Mix.displayFit(timeDomainSignal,t)
                
                % If user wants to keep component, add it to the object
                output_str1 = lower(input('Keep fit? [y/n]','s'));
                keepFit = length(findstr(output_str1, 'y')>0);
                if(keepFit)
                    obj.amp = temp_NMR_Mix.amp;
                    obj.freq = temp_NMR_Mix.freq;
                    obj.fwhm = temp_NMR_Mix.fwhm;
                    obj.phase = temp_NMR_Mix.phase;
                end
                
                % If user wants to keep component, add it to the object
                output_str2 = lower(input('\tFit more peaks? [y/n]','s'));
                continueFitting = length(findstr(output_str2, 'y')>0);
            end
        end
        
        function obj = fitTimeDomainSignal(obj, timeDomainSignal, t)
            % Fits exponentially decaying components to the given time
            % domain signal
            
            % Convert to the matrix form
            guess_compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            % Fit time domain signal
            fitfunc = @(compMat,t)sum(obj.calcTimeDomainSignalFromMatrix(compMat,t),2);
            ub = []; % Lower bounds
            lb = []; % Upper bounds
            fitoptions = optimoptions('lsqcurvefit');
            %                         fitoptions.Display = 'iter-detailed';
            fitoptions.MaxIter = 400;
            fitoptions.TolFun=1E-900;
            fitoptions.TolX = 1E-20;
            fitoptions.MaxFunEvals = 500;
            [fit_compMat,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqcurvefit(fitfunc,guess_compMat,t,timeDomainSignal,lb,ub,fitoptions);
            
            % Convert back to physical units and save to object
            [obj.amp, obj.freq, obj.fwhm, obj.phase] = ...
                obj.convertMatrixToPhysical(fit_compMat);
        end
        
        function displayFit(obj, timeDomainSignal,t)
            % Calculate spectrum
            dwell_time = t(2)-t(1);
            f = NMR_Mix.calcFftFreq(t);
            measuredSpectrum = dwell_time*fftshift(fft(timeDomainSignal));
            
            % Calculate fitted and residual spectrums
            fittedSpectrum = sum(obj.calcSpectralDomainSignal(f),2);
            residualSpectrum = measuredSpectrum - fittedSpectrum;
            
            % Calculate lorentzian curves for each component
            lorentzianCurves = obj.calcLorentzianCurves(f);
            nComponents = length(obj.amp);
            fMat = repmat(f,[1 nComponents]);
            
            legendStrings = cell(1, nComponents);
            for iComp=1:nComponents
                legendStrings{iComp} = ['C\_' sprintf('%03.0f',iComp)];
            end
            
            % Show results to user
            ax2 = subplot(4,1,2);
            plot(f,abs(measuredSpectrum),'-b');
            hold on;
            plot(f,abs(fittedSpectrum),'-g');
            plot(f,abs(residualSpectrum),'-r');
            hold off;
            ylabel('Magnitude Intensity');
            set(ax2,'xticklabel',{[]}) ;
            set(ax2,'XDir','reverse');
            
            ax3 = subplot(4,1,3);
            plot(f,real(measuredSpectrum),'-b');
            hold on;
            plot(f,real(fittedSpectrum),'-g');
            plot(f,real(residualSpectrum),'-r');
            hold off;
            ylabel('Real Intensity');
            set(ax3,'xticklabel',{[]});
            set(ax3,'XDir','reverse');
            
            ax4 = subplot(4,1,4);
            plot(f,imag(measuredSpectrum),'-b');
            hold on;
            plot(f,imag(fittedSpectrum),'-g');
            plot(f,imag(residualSpectrum),'-r');
            hold off;
            xlabel('Spectral Frequency (Hz)');
            ylabel('Imaginary Intensity');
            legend('Measured','Fitted','Residual');
            set(ax4,'XDir','reverse');
            
            if(~isempty(obj.fref))
                % Add PPM axis
                ax1ppm = subplot(4,1,1);
                set(ax1ppm,'units','normalized',...
                    'XAxisLocation','top','YAxisLocation','right',...
                    'YTick',[],'YTickLabel',[],'Color','none');
                
                ax1 = axes();
            else
                ax1 = subplot(4,1,1);
            end
                        
            hold off;
            plot(fMat,lorentzianCurves);
            legend(legendStrings);
            ylabel('Component Intensity');
            set(ax1,'xticklabel',{[]}) ;
            set(ax1,'XDir','reverse');
            if(~isempty(obj.fref))
                set(ax1ppm,'XDir','reverse');
            end
            
            % Clean up gaps between plots
            if(~isempty(obj.fref))
                ax1Pos = get(ax1ppm,'Position');
            else
                ax1Pos = get(ax1,'Position');
            end
            ax2Pos = get(ax2,'Position');
            ax3Pos = get(ax3,'Position');
            ax4Pos = get(ax4,'Position');
            
            % Better position each axes object
            set(ax1,'Position',[ax1Pos(1) ax1Pos(2) ax1Pos(3) 0.2]);
            set(ax2,'Position',[ax2Pos(1) ax2Pos(2) ax2Pos(3) 0.2]);
            set(ax3,'Position',[ax3Pos(1) ax3Pos(2) ax3Pos(3) 0.2]);
            set(ax4,'Position',[ax4Pos(1) ax4Pos(2) ax4Pos(3) 0.2]);
            
            % Keep all x axes in sinc
            linkaxes([ax1,ax2,ax3,ax4],'x');
            
            if(~isempty(obj.fref))
                % Set size of ppm axis
                set(ax1ppm,'Position',[ax1Pos(1) ax1Pos(2) ax1Pos(3) 0.2]);
                
                % Initialize XLim in correct units
                set(ax1ppm,'xlim',NMR_Mix.relFreqToPpm(get(ax2,'XLim'),obj.fref));
                
                % Keep ppm axiz in sinc with freq axis
                xLimListener = addlistener( ax1, 'XLim', 'PostSet', ...
                    @(src,evt) set(ax1ppm,'XLim',...
                    NMR_Mix.relFreqToPpm(get(ax1,'XLim'),obj.fref)) );
            end
        end
    end
    
    % These are convenience methods for unit conversions, etc
    methods (Static=true)
        function w = freqToAngFreq(f)
            w = 2*pi*f;
        end
        function f = AngFreqToFreq(w)
            f = w./(2*pi);
        end
        function ppm = relFreqToPpm(f,fref)
            % Calculates PPM values from a relative frequency (difference
            % from reference) and reference frequency.
            ppm = 1E6*(f./fref);
        end
        function ppm = freqToPpm(f,fref)
            % Calculates PPM values from a frequency and reference
            % frequency.
            ppm = 1E6*((f - fref)./fref);
        end
        function f = ppmToFreq(ppm,fref)
            % Calculates frequencies from PPM values and reference
            % frequency
            f = ppm*1E-6*fref+fref;
        end
        function f = ppmToRelfreq(ppm,fref)
            % Calculates relative frequencies (difference from referecne
            % frequency) from PPM values and reference frequency
            f = ppm*1E-6*fref;
        end
        function f = calcFftFreq(t)
            nSamples = length(t);
            dwell_time = t(2)-t(1);
            f = linspace(-0.5,0.5,nSamples+1)/dwell_time;
            f = f(1:(end-1))'; % Take off last sample to have nSamples
        end
        function t = calcIfftTime(f)
            nSamples = length(f);
            freq_gap = f(2)-f(1);
            t = linspace(-0.5,0.5,nSamples+1)/freq_gap;
            t = t(1:(end-1)); % Take off last sample to have nSamples
        end
    end
    
    methods (Access=protected)
        function componentMatrix = convertPhysicalToMatrix(obj, ...
                amplitude, frequency, fwhm, phase)
            % Converts from physical parameters to the more compact complex
            % notation. This is useful for fitting purposes.
            time_invariant = amplitude.*exp(1i*degtorad(phase)); % time invariant component
            time_variant = (-pi*fwhm + 1i*obj.freqToAngFreq(frequency)); % time variant component
            componentMatrix = [time_invariant time_variant];
        end
        
        function [amplitude, frequency, fwhm, phase] = ...
                convertMatrixToPhysical(obj, componentMatrix)
            % Converts from the compact complex notation to physical units.
            amplitude   = abs(componentMatrix(:,1));        % arb units
            phase = radtodeg(angle(componentMatrix(:,1)));		% Degrees
            fwhm  = real(componentMatrix(:,2))/-pi;			% Hz
            frequency  = obj.AngFreqToFreq(imag(componentMatrix(:,2))); % Hz
        end
        
        function timeDomainSignal = ...
                calcTimeDomainSignalFromMatrix(obj,componentMatrix,t)
            % Calculates the time domain signal from the mix of NMR
            % components in the component matrix
            nComponents = size(componentMatrix,1);
            timeDomainSignal = zeros(length(t),nComponents);
            for iComp = 1:nComponents
                timeDomainSignal(:,iComp) =	componentMatrix(iComp,1) .* ...
                    exp(componentMatrix(iComp,2).*t);
            end
        end
        
        function spectralDomainSignal = ...
                calcSpectralDomainSignalFromMatrix(obj,componentMatrix,f)
            % Calculates the spectral domain signal from the mix of NMR
            % components in the component matrix
            nComponents = size(componentMatrix,1);
            spectralDomainSignal = zeros(length(f),nComponents);
            for iComp = 1:nComponents
                spectralDomainSignal(:,iComp) =	...
                    componentMatrix(iComp,1) ./ ...
                    (1i*obj.freqToAngFreq(f)-componentMatrix(iComp,2));
            end
        end
        
        function lorentzianCurves = calcLorentzianCurvesFromMatrix(obj, componentMatrix, f)
            % Calculates unphased lorentzians from the mix of NMR
            % components in the component matrix
            nComponents = size(componentMatrix,1);
            lorentzianCurves = zeros(length(f),nComponents);
            for iComp = 1:nComponents
                lorentzianCurves(:,iComp) = ...
                    real(abs(componentMatrix(iComp,1)) ./ ...
                    (1i*obj.freqToAngFreq(f)-componentMatrix(iComp,2)));
            end
        end
        
        function dispersiveCurves = calcDispersiveCurvesFromMatrix(obj, componentMatrix, f)
            % Calculates unphased lorentzians from the mix of NMR
            % components in the component matrix
            nComponents = size(componentMatrix,1);
            dispersiveCurves = zeros(length(f),nComponents);
            for iComp = 1:nComponents
                dispersiveCurves(:,iComp) = ...
                    imag(abs(componentMatrix(iComp,1)) ./ ...
                    (1i*obj.freqToAngFreq(f)-componentMatrix(iComp,2)));
            end
        end
        
        function header = getHeader(obj)
            nComponents = length(obj.amp);
            if(~isempty(obj.fref))
                header = ['Component:  |  Amp(arb) | Freq(Hz) | Freq(ppm) | FWHM(Hz) | Phase(deg) |'];
            else
                header = ['Component:  |  Amp(arb) | Freq(Hz) |  FWHM(Hz) | Phase(deg) |'];
            end
            
        end
        
        function propgrp = getPropertyGroups(obj)
            nComp = length(obj.amp);
            propList = struct();
            for iComp=1:nComp
                name = ['C_' sprintf('%03.0f',iComp)];
                details = [ '| ' sprintf('%8.3e',obj.amp(iComp)) ' | ' ...
                    sprintf('%+8.2f',obj.freq(iComp)) ' | '];
                if(~isempty(obj.fref))
                    details = [details sprintf('%+9.2f',NMR_Mix.relFreqToPpm(obj.freq(iComp),obj.fref)) ' | '];
                end
                details = [details sprintf('%8.2f',obj.fwhm(iComp)) ' | ' ...
                    sprintf('%+9.2f',obj.phase(iComp))];
                propList = setfield(propList,name,details);
            end
            propgrp = matlab.mixin.util.PropertyGroup(propList);
        end
    end
end
