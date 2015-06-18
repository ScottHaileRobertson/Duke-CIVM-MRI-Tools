%% NMR_Fit
% A class that fits a series of exponentially decaying components to a
% NMR_Mix object. The class knows how to fit to a time domain signal with
% or without constraints.
%
classdef NMR_Fit
    
    properties
        nmrMix;
        ub;
        lb;
        t;
        f;
        timeDomainSignal;
        spectralDomainSignal;
        zeroPadding;
        lineBroadening;
    end
    
    methods
        function obj = NMR_Fit(time_domain_signal, t, zero_padding, ...
                line_broadening, amp, freq, fwhm, phase)
            obj.t = t(:);
            nSamples = length(obj.t);
            obj.zeroPadding = zero_padding;
            if(isempty(obj.zeroPadding))
                obj.zeroPadding = nSamples;
            end
            obj.lineBroadening = line_broadening;
            if(isempty(obj.lineBroadening))
                obj.lineBroadening = 0;
            end
            
            % Apply linebroadening
            obj.timeDomainSignal = time_domain_signal.*exp(-obj.lineBroadening*obj.t);
            
            % Calculate spectrum
            dwell_time = (obj.t(2)-obj.t(1));
            obj.spectralDomainSignal = dwell_time...
                *fftshift(fft(obj.timeDomainSignal,obj.zeroPadding));
            
            % Calculate freq samples
            obj.f = linspace(-0.5,0.5,obj.zeroPadding+1)/dwell_time;
            obj.f = obj.f(1:(end-1))'; % Take off last sample to have nSamples
            
            % Create a NMR_Mix object
            obj.nmrMix = NMR_Mix(amp,freq,fwhm,phase);
        end
        
        function obj = setBounds(obj, amp_lb, amp_ub, freq_lb, freq_ub,...
                fwhm_lb, fwhm_ub, phase_lb, phase_ub)
            % Sets the constraints for fitting this NMR_Mix object. Be very
            % careful to only add bounds just before you fit the signal.
            % Adding additional components to the nmr_mix will erase the
            % bounds, as they will potentially get out of order.
            obj.ub = [amp_ub(:)'; freq_ub(:)'; fwhm_ub(:)'; phase_ub(:)'];
            obj.lb = [amp_lb(:)'; freq_lb(:)'; fwhm_lb(:)'; phase_lb(:)'];
            
            % Check that constraints are possible
            if(any(obj.ub < obj.lb))
                error('Impossible constraint set.');
            end
        end
        
        function [fit_amp, fit_freq, fit_fwhm, fit_phase] = ...
                fitComponentToResidual(obj)
            % Fits a single exponentially decaying component to the
            % residual signal. This is useful in peak fitting. Note,
            % however, that this function only fits based on residuals,
            % thus you are fitting the residuals, and not the signal. After
            % fitting the residuals, it is recommended that you "refit" all
            % the components as they will likely change slightly to better
            % accomodate the new compnent.
            
            % If there are already components, fit them first
            if(~isempty(obj.nmrMix.amp))
                currentSpectrum = obj.spectralDomainSignal - ...
                    obj.nmrMix.calcSpectralDomainSignal(obj.t);
            else
                currentSpectrum = obj.spectralDomainSignal;
            end
            
            % Find the frequency of the next peak
            [maxVal maxIdx] = max(abs(currentSpectrum));
            peakFreq = obj.f(maxIdx);
            
            % Use time domain curve fitting to fit residual signal
            residualTimeDomainSignal = obj.timeDomainSignal - obj.nmrMix.calcTimeDomainSignal(obj.t);
            largestPeakComponent = NMR_Fit(residualTimeDomainSignal, obj.t, ...
                obj.zeroPadding,obj.lineBroadening, 1, peakFreq, 100, 0);
            [fit_amp, fit_freq, fit_fwhm, fit_phase] = ...
                largestPeakComponent.calcTimeDomainSignalFit();
        end
        
        function obj = autoAddComponent(obj)
            % This function attemps to automatically fit and add a new
            % exponentially decaying signal to the NMR_mix object
            
            % Fit next component to residual signal, then add it to this
            % NMR_Mix object
            [add_amp, add_freq, add_fwhm, add_phase] = ...
                obj.fitComponentToResidual();
            
            % Add fitted component to NMR_mix
            obj.nmrMix = obj.nmrMix.addComponents(add_amp, add_freq, add_fwhm, add_phase);
            
            % Refit all components after addition of latest component
            obj = obj.fitTimeDomainSignal();
            
            if(~isempty(obj.ub))
                % Remove bounds in case sorted order changed
                obj.ub = [];
                obj.lb = [];
                warning('Deleting bounds - you need to call setBounds again with the new components bounds');
            end
        end
        
        function obj = autoAddComponents(obj, nComponents)
            % This function attemps to automatically fit and add
            % n (defined by nComponents) new exponentially decaying signals
            % to the NMR_mix object
            
            for iComp = 1:nComponents
                obj = obj.autoAddComponent();
            end
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
                % Create temporary NRM_Fix object and add newest component
                tempFit = obj;
                tempFit = tempFit.autoAddComponent();
                
                % Show fit
                tempFit.plotFit();
                
                % Report new fits
                tempFit.describe();
                
                % If user wants to keep component, add it to the object
                output_str1 = lower(input('Keep fit? [y/n]','s'));
                keepFit = length(findstr(output_str1, 'y')>0);
                if(keepFit)
                    obj = tempFit;
                end
                
                % If user wants to fit more peaks, continue fitting...
                output_str2 = lower(input('\tFit more peaks? [y/n]','s'));
                continueFitting = length(findstr(output_str2, 'y')>0);
            end
        end
        
        function [fit_amp, fit_freq, fit_fwhm, fit_phase] = calcTimeDomainSignalFit(obj)
            % Fits exponentially decaying components to the given time
            % domain signal and provides the results without saving them to
            % this object. To save the results to the NMR_Mix, use
            % fitTimeDomainSignal
            
            % Put all components into a single matrix
            guess = [obj.nmrMix.amp; obj.nmrMix.freq; ...
                obj.nmrMix.fwhm; obj.nmrMix.phase];
            
            fitoptions = optimoptions('lsqcurvefit');
            %             fitoptions.Display = 'iter-detailed';
            %             fitoptions.Display = 'final-detailed';
            fitoptions.Display = 'off';
            fitoptions.MaxIter = 1000;
            fitoptions.TolFun=1E-900;
            fitoptions.TolX = 1E-8;
            fitoptions.FinDiffType = 'central';
            fitoptions.Algorithm = 'trust-region-reflective';
            fitoptions.MaxFunEvals = 5000;
            fit_params = lsqcurvefit(@obj.calcTimeSigRealImag,guess,obj.t(:),...
                [real(obj.timeDomainSignal(:)),imag(obj.timeDomainSignal(:))],...
                obj.lb(:),obj.ub(:),fitoptions);
            
            % Separate out the components from the matrix
            fit_vec = fit_params(1,:).*exp(1i*pi*fit_params(4,:)/180);
            fit_amp = abs(fit_vec);
            fit_freq = fit_params(2,:);
            fit_fwhm = fit_params(3,:);
            fit_phase = angle(fit_vec)*180/pi;
        end
        
        function obj = fitTimeDomainSignal(obj)
            % Fits exponentially decaying components to the given time
            % domain signal and saves the results to this NMR_Mix object.
            % To just return the fits and not save the results, use
            % calcTimeDomainSignalFit
            [fit_amp, fit_freq, fit_fwhm, fit_phase] = obj.calcTimeDomainSignalFit();
            
            % Save fits
            obj.nmrMix = obj.nmrMix.resetComponents(fit_amp, fit_freq, fit_fwhm, fit_phase);
        end
        
        function realImagSig = calcTimeSigRealImag(obj,nmr_params,t)
            % A function used in fitting to allow constraints for complex
            % fitting. This is the same as calling calcTimeDomainSignal of
            % the NMR_Mix object, only it separates the real and imaginary
            % parts to allow for constraints to be used in fitting.       
              tmpNmrMix = NMR_Mix(nmr_params(1,:), nmr_params(2,:), ...
                  nmr_params(3,:), nmr_params(4,:));
              complexSig = sum(tmpNmrMix.calcTimeDomainSignal(t),2);
              realImagSig = [real(complexSig) imag(complexSig)];
        end
        
        function plotFit(obj)
            % Calculate fitted and residual spectrums
            individualSpectrums = obj.nmrMix.calcComponentSpectralDomainSignal(obj.f);
            fittedSpectrum = obj.nmrMix.calcSpectralDomainSignal(obj.f);
            residualSpectrum = obj.spectralDomainSignal - fittedSpectrum;
            
            % Calculate lorentzian curves for each component
            nComponents = length(obj.nmrMix.amp);
            fMat = repmat(obj.f,[1 nComponents]);
            
            legendStrings = cell(1, nComponents);
            for iComp=1:nComponents
                legendStrings{iComp} = ['C\_' sprintf('%03.0f',iComp)];
            end
            
            % Show results to user
            ax2 = subplot(5,1,2);
            plot(obj.f,abs(obj.spectralDomainSignal),'-b');
            hold on;
            plot(obj.f,abs(fittedSpectrum),'-g');
            plot(obj.f,abs(residualSpectrum),'-r');
            hold off;
            ylabel('Magnitude Intensity');
            set(ax2,'xticklabel',{[]}) ;
            set(ax2,'XDir','reverse');
            
            ax3 = subplot(5,1,3);
            phaseSig = angle(obj.spectralDomainSignal);
            phaseFit = angle(fittedSpectrum);
            phaseDelta = phaseSig - phaseFit;
            plot(obj.f, phaseSig,'-b');
            hold on;
            plot(obj.f,phaseFit,'-g');
            plot(obj.f,phaseDelta,'-r');
            
            hold off;
            ylabel('Phase (Radians)');
            set(ax3,'xticklabel',{[]}) ;
            set(ax3,'XDir','reverse');
            
            ax4 = subplot(5,1,4);
            plot(obj.f,real(obj.spectralDomainSignal),'-b');
            hold on;
            plot(obj.f,real(fittedSpectrum),'-g');
            plot(obj.f,real(residualSpectrum),'-r');
            hold off;
            ylabel('Real Intensity');
            set(ax4,'xticklabel',{[]});
            set(ax4,'XDir','reverse');
            
            ax5 = subplot(5,1,5);
            plot(obj.f,imag(obj.spectralDomainSignal),'-b');
            hold on;
            plot(obj.f,imag(fittedSpectrum),'-g');
            plot(obj.f,imag(residualSpectrum),'-r');
            hold off;
            xlabel('Spectral Frequency (Hz)');
            ylabel('Imaginary Intensity');
            legend('Measured','Fitted','Residual');
            set(ax5,'XDir','reverse');
            
            %             if(~isempty(obj.fref))
            %                 % Add PPM axis
            %                 ax1ppm = subplot(5,1,1);
            %                 set(ax1ppm,'units','normalized',...
            %                     'XAxisLocation','top','YAxisLocation','right',...
            %                     'YTick',[],'YTickLabel',[],'Color','none');
            %
            %                 ax1 = axes('Position',get(ax1ppm,'Position'));
            %             else
            ax1 = subplot(5,1,1);
            %             end
            
            plot(fMat,real(individualSpectrums));
            legend(legendStrings);
            ylabel('Component Intensity');
            set(ax1,'xticklabel',{[]}) ;
            set(ax1,'XDir','reverse');
%             if(~isempty(obj.fref))
%                 set(ax1ppm,'XDir','reverse');
%             end
            
            % Keep all x axes in sinc
            linkaxes([ax1,ax2,ax3,ax4 ax5],'x');
            
            %             if(~isempty(obj.fref))
            %                 % Initialize XLim in correct units
            %                 set(ax1ppm,'xlim',NMR_Mix.relFreqToPpm(get(ax2,'XLim'),obj.fref));
            %
            %                 % Keep ppm axiz in sinc with freq axis
            %                 xLimListener = addlistener( ax1, 'XLim', 'PostSet', ...
            %                     @(src,evt) set(ax1ppm,'XLim',...
            %                     NMR_Mix.relFreqToPpm(get(ax1,'XLim'),obj.fref)) );
            %             end
        end
        
        function describe(obj)
            disp('Amp (arbs)  Freq (Hz)  Linewidth(Hz)  Phase(degrees)');
            for iComp = 1:length(obj.nmrMix.amp)
                disp([sprintf('%8.3e',obj.nmrMix.amp(iComp)) ' ' ...
                    sprintf('%+8.2f',obj.nmrMix.freq(iComp))  '  ' ...
                    sprintf('%8.2f',obj.nmrMix.fwhm(iComp)) '  ' ...
                    sprintf('%+9.2f',obj.nmrMix.phase(iComp))]);
            end
        end
    end
end