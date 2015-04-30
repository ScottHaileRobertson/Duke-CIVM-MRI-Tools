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
        timeDomainSignal;
        t;
        spectralDomainSignal; 
        f; 
        
        zeroPadSize;
        linebroadening;
        
        amp;   % Amplitude of each component (arbs)
        freq;  % Frequency of each component in Hz/Rad
        fwhm;  % FWHM of each component in Hz
        phase; % Phase of each component in radians
        fref = [];
    end
    
    methods
        function obj = NMR_Mix(time_domain_signal, time_vec, varargin)
            % Constructs an NMR_Mix object with the given time domain
            % signal  and optional knowledge of component amplitudes
            % (arbs), frequencies (in Hz), fwhms (in Hz), and phases (in
            % degrees). To start with an empty object (ex for fitting),
            % just use NMR_mix(time_domain_signal, time_vec)
            obj.t = time_vec(:);
            obj.timeDomainSignal = time_domain_signal(:);
            
            if(nargin > 2)
                obj.amp = varargin{1}(:);
                obj.freq = varargin{2}(:); % In Hz
                obj.fwhm = varargin{3}(:);       % In Hz
                obj.phase = varargin{4}(:);     % In Degrees
                
                % Resort object by frequencies
                obj = obj.sortByFreq();
                if(nargin > 6)
                    % Populate zeropadding(default is no padding)
                    obj.zeroPadSize = varargin{5};
                    
                    % Populate line broadening (default is none)
                    obj.linebroadening = varargin{6};
                    
                    if(nargin == 9)
                        obj.fref = varargin{7};
                    end
                end
                if(isempty(obj.zeroPadSize))
                    obj.zeroPadSize = length(obj.t);
                end
                if(isempty(obj.linebroadening))
                    obj.linebroadening = 0;
                end
                
                % Calculate spectrum\
                nSamples = length(obj.t);
                dwell_time = (obj.t(2)-obj.t(1));
                obj.spectralDomainSignal = dwell_time...
                    *fftshift(fft(obj.timeDomainSignal,obj.zeroPadSize));
                
                % Calculate Frequencies
                obj.f = linspace(-0.5,0.5,obj.zeroPadSize+1)/dwell_time;
                obj.f = obj.f(1:(end-1))'; % Take off last sample to have nSamples
            end
            
            % Apply linebroadening
            obj.timeDomainSignal = obj.timeDomainSignal.*exp(-obj.linebroadening*obj.t);
        end
        
        function obj = addComponent(obj,amplitude, frequency, fwhm, phase)
            obj.amp = [obj.amp; amplitude];
            obj.freq = [obj.freq; frequency];
            obj.fwhm = [obj.fwhm; fwhm];
            obj.phase = [obj.phase; phase];
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = removeComponents(obj, componentNumbers)
            obj.amp(componentNumbers) = [];
            obj.freq(componentNumbers) = [];
            obj.fwhm(componentNumbers) = [];
            obj.phase(componentNumbers) = [];
        end
        
        function obj = sortByFreq(obj)
            % Sort frequencies
            [sortFreq sortIdx] = sort(obj.freq, 'descend');
            
            % resort object
            obj.amp = obj.amp(sortIdx);
            obj.freq = obj.freq(sortIdx);
            obj.fwhm = obj.fwhm(sortIdx);
            obj.phase = obj.phase(sortIdx);
        end
        
        function timeDomainSignal = calcTimeDomainSignal(obj,varargin)
            % Calculates the time domain signal from the mix of NMR
            % components at the given time points (t is in sec). Note, this
            % function returns the time signal for each individual
            % component. The overall "mix" signal can be obtained by
            % summing accross all components.
            
            % Convert to the matrix form
            compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            % Get time vector from user or object if not given
            if(nargin == 2)
                t = varargin{1}(:);
            else
                t = obj.t;
            end
            
            % Calculate time domain signal using matrix form
            timeDomainSignal= obj.calcTimeDomainSignalFromMatrix(compMat,t);
        end
        
        function spectralDomainSignal = calcSpectralDomainSignal(obj,varargin)
            % Calculates the spectral domain signal from the mix of NMR
            % components at the given spectral frequencies (f is in Hz).
            % Note, this function returns the spectrum for each individual
            % component. The overall "mix" spectrum can be obtained by
            % summing accross all components.
            
            % Convert to the matrix form
            compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            % Get frequency vector from user or object if not given
            if(nargin == 2)
                f = varargin{1}(:);
            else
                f = obj.f;
            end
            
            % Calculate spectral domain signal using matrix form
            spectralDomainSignal = obj.calcSpectralDomainSignalFromMatrix(compMat,f);
        end
        
        function lorentzianCurves = calcLorentzianCurves(obj,varargin)
            % Calculates unphased lorenzian curves from the mix of NMR
            % components at the given spectral frequencies (f is in Hz).
            % Note, this function returns a curve for each individual
            % component. The overall "mix" of lorentzian curves can be
            % obtained by summing accross all components.
            
            % convert to matrix form
            compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            % Get frequency vector from user or object if not given
            if(nargin == 2)
                f = varargin{1}(:);
            else
                f = obj.f;
            end
            
            lorentzianCurves = obj.calcLorentzianCurvesFromMatrix(compMat, f);
        end
        
        function dispersiveCurves = calcDispersiveCurves(obj,varargin)
            % Calculates unphased dispersive curves from the mix of NMR
            % components at the given spectral frequencies (f is in Hz).
            % Note, this function returns a curve for each individual
            % component. The overall "mix" of dispersive curves can be
            % obtained by summing accross all components.
            
            % convert to matrix form
            compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
           
            % Get frequency vector from user or object if not given
            if(nargin == 2)
                f = varargin{1}(:);
            else
                f = obj.f;
            end
            
            dispersiveCurves = obj.calcDispersiveCurvesFromMatrix(compMat, f);
        end
        
        function [amp freq fwhm phase] = fitComponentToResidual(obj)
            % Fits a single exponentially decaying component to the
            % residual signal. This is useful in peak fitting. Note,
            % however, that this function only fits based on residuals,
            % thus you are fitting the residuals, and not the signal. After
            % fitting the residuals, it is recommended that you "refit" all
            % the components as they will likely change slightly to better
            % accomodate the new compnent.
            
            % If there are already components, fit them first
            if(~isempty(obj.amp) | ~isempty(obj.freq) | ...
                    ~isempty(obj.fwhm) | ~isempty(obj.phase))
                currentSpectrum = obj.spectralDomainSignal - ...
                    sum(obj.calcSpectralDomainSignal(),2);
            else
                currentSpectrum = obj.spectralDomainSignal;
            end
            
            % Find the frequency of the next peak
            [maxVal maxIdx] = max(abs(currentSpectrum));
            peakFreq = obj.f(maxIdx);
            
            % Find the phase of this peak
            nPhases = 360;
            phases = linspace(-pi,pi,nPhases+1);
            phases = phases(1:(end-1));
            [maxVal maxIdx] = max(real(currentSpectrum(maxIdx).*exp(1i*phases)));
            peakPhase = 180*phases(maxIdx)/pi; % in deg
            
%             phaseMat = repmat(phases,[length(currentSpectrum) 1]);
%             specMat = repmat(currentSpectrum,[1 nPhases]);
%             specMat = specMat.*phaseMat;
            
            % Use curve fitting to fit residual signal
            residualTimeDomainSignal = obj.timeDomainSignal - sum(obj.calcTimeDomainSignal(obj.t),2);
            largestPeakComponent = NMR_Mix(residualTimeDomainSignal,obj.t,...
                1,peakFreq,10,peakPhase,obj.zeroPadSize,obj.linebroadening,obj.fref); %Guess amplitude and FWHM for now
            largestPeakComponent = largestPeakComponent.fitTimeDomainSignal();
            
            % Return fit parameters
            amp   = largestPeakComponent.amp;
            freq  = largestPeakComponent.freq;
            fwhm  = largestPeakComponent.fwhm;
            phase = largestPeakComponent.phase;
        end
        
        function obj = autoAddComponent(obj)
            % This function attemps to automatically fit and add a new
            % exponentially decaying signal to the NMR_mix object
            
            % Fit next component to residual signal, then add it to this
            % NMR_Mix object
            [add_amp add_freq, add_fwhm add_phase] = ...
                obj.fitComponentToResidual();
            
            % Add fitted component to NMR_mix
            obj.amp = [obj.amp; add_amp];
            obj.freq = [obj.freq; add_freq];
            obj.fwhm = [obj.fwhm; add_fwhm];
            obj.phase = [obj.phase; add_phase];
            
            % Refit all components after addition of latest component
            obj = obj.fitTimeDomainSignal();
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = autoAddComponents(obj, nComponents)
            % This function attemps to automatically fit and add
            % n (defined by nComponents) new exponentially decaying signals
            % to the NMR_mix object
            
            nIter = nComponents - length(obj.amp);
            for iComp = 1:nIter
                obj = obj.autoAddComponent();
            end
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = fitTool(obj)
            % This function attemps to adds exponentially decaying signal
            % components until the user stops the fit.
                        
            % Create a figure to show results
            thisFigure = gcf();
            clf;
            
            % Keep fitting until user says stop
            continueFitting = true;
            while(continueFitting)
                % Fit next component to residual signal, then add it to this
                % NMR_Mix object
                [add_amp add_freq, add_fwhm add_phase] = obj.fitComponentToResidual();
                
                % Create temporary NRM_mix object and add newest component
                temp_NMR_Mix = NMR_Mix(obj.timeDomainSignal, obj.t,...
                    obj.amp, obj.freq, obj.fwhm, obj.phase, ...
                    obj.zeroPadSize, obj.linebroadening, obj.fref);
                temp_NMR_Mix.amp = [temp_NMR_Mix.amp; add_amp];
                temp_NMR_Mix.freq = [temp_NMR_Mix.freq; add_freq];
                temp_NMR_Mix.fwhm = [temp_NMR_Mix.fwhm; add_fwhm];
                temp_NMR_Mix.phase = [temp_NMR_Mix.phase; add_phase];
                
                % Refit all components after addition of latest component
                temp_NMR_Mix = temp_NMR_Mix.fitTimeDomainSignal();
                
                % Show fit
                temp_NMR_Mix.plotFit();
                
                % Resort object by frequencies
                obj = obj.sortByFreq();
                
                % Report new fits
                obj
                
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
        
        function obj = fitTimeDomainSignal(obj)
            % Fits exponentially decaying components to the given time
            % domain signal
            
            % Convert to the matrix form
            guess_compMat = obj.convertPhysicalToMatrix(obj.amp, ...
                obj.freq, obj.fwhm, obj.phase);
            
            % Fit time domain signal
            fitfunc = @(input_params,t)sum(obj.calcTimeDomainSignalFromMatrix(input_params,t),2);
            
            ub = []; % Lower bounds
            lb = []; % Upper bounds
            fitoptions = optimoptions('lsqcurvefit','Display','off');
            %                         fitoptions.Display = 'iter-detailed';
%             fitoptions.MaxIter = 5000;
%             fitoptions.TolFun=1E-900;
%             fitoptions.TolX = 1E-900;
%             fitoptions.FinDiffType = 'central';
%             fitoptions.MaxFunEvals = 1000;
            [fit_params,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqcurvefit(fitfunc,guess_compMat,obj.t,...
                obj.timeDomainSignal,lb,ub,fitoptions);
            
            % Convert back to physical units and save to object
            [obj.amp, obj.freq, obj.fwhm, obj.phase] = ...
                obj.convertMatrixToPhysical(fit_params);
        end
        
        function plotFit(obj)
            % Calculate fitted and residual spectrums
            fittedSpectrum = sum(obj.calcSpectralDomainSignal,2);
            residualSpectrum = obj.spectralDomainSignal - fittedSpectrum;
            
            % Calculate lorentzian curves for each component
            lorentzianCurves = obj.calcLorentzianCurves();
            nComponents = length(obj.amp);
            fMat = repmat(obj.f,[1 nComponents]);
            
            legendStrings = cell(1, nComponents);
            for iComp=1:nComponents
                legendStrings{iComp} = ['C\_' sprintf('%03.0f',iComp)];
            end
            
            % Show results to user
            ax2 = subplot(4,1,2);
            plot(obj.f,abs(obj.spectralDomainSignal),'-b');
            hold on;
            plot(obj.f,abs(fittedSpectrum),'-g');
            plot(obj.f,abs(residualSpectrum),'-r');
            hold off;
            ylabel('Magnitude Intensity');
            set(ax2,'xticklabel',{[]}) ;
            set(ax2,'XDir','reverse');
            
            ax3 = subplot(4,1,3);
            plot(obj.f,real(obj.spectralDomainSignal),'-b');
            hold on;
            plot(obj.f,real(fittedSpectrum),'-g');
            plot(obj.f,real(residualSpectrum),'-r');
            hold off;
            ylabel('Real Intensity');
            set(ax3,'xticklabel',{[]});
            set(ax3,'XDir','reverse');
            
            ax4 = subplot(4,1,4);
            plot(obj.f,imag(obj.spectralDomainSignal),'-b');
            hold on;
            plot(obj.f,imag(fittedSpectrum),'-g');
            plot(obj.f,imag(residualSpectrum),'-r');
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
            frequency  = obj.AngFreqToFreq(imag(componentMatrix(:,2))); % Hz
            fwhm  = real(componentMatrix(:,2))/-pi;			% Hz
            phase = radtodeg(angle(componentMatrix(:,1)));		% Degrees
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
