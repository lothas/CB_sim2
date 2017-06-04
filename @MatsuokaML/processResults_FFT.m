function [periodsFFT,harmonics_amp,harmonics_phase, signal_amp] = processResults_FFT( obj, CPG_out,...
    T, graphicsFlag)
    % this function is for the 2N symmetric Matsuoka CPG case only
    
    % INPUTS:
    % 'T' - time vector
    % 'CPG_out' - the output of the CPG. defined as 'x2-x1'
    
    % OUTPUTS:
    % *) 'periodsFFT' - the periods of the signal (==1/harmonicsFreq)
    % *) 'harmonics_amp' - the amplitude of each harmonics
    % *) 'harmonics_phase' - the phase of each harmonics
    % *) 'signal_amp' - the amplitude of the CPG output
    
    % process parameters:
    percent_of_sig2process = 0.25; % aliminate the transient part
    
    X = CPG_out;
    T = T';
    
    % Turn off findpeaks warning
    warning('off','signal:findpeaks:largeMinPeakHeight');

    dt = mean(diff(T));
    fs = 1/dt;
    
    signal_start = floor(percent_of_sig2process * length(T));
    longest_signal = 1; signal_length = 1;
    
    signal2Process = X(1,signal_start:end);
    % aliminate the 'zero' frequency (bias)
    signal2Process = signal2Process - min(signal2Process);
    
    % check if the signal is non periodic
    if std(signal2Process) < 1e-3 % check standart deviation:
        fundanetalHarmonics = NaN;
        % TODO: implement amplitude detection
        signal_amp = []; % mean(signals'); %#ok<UDIM>
        return;
    end

    % find the peaks in the time domain in order to
    %   take a 'N' number of periods. 
    [~,sigPeak_locs]=findpeaks(signal2Process, ...
                'MinPeakheight',...
                0.85*max(signal2Process));
    
    % check if we get more than one peak        
    if length(sigPeak_locs) > 1
        peak_locs(1,:) = [sigPeak_locs(1), sigPeak_locs(end)];
    else
        fundanetalHarmonics = NaN;
        signal_amp = mean(X'); %#ok<UDIM>
        return
    end
    % get the indices of the original signal:
    locs = signal_start+[peak_locs(1,1), ...
        peak_locs(1,2)];
    
    % prepare to perform the FFT:
    signal2FFT = X(:,locs(1):locs(end));

    m = length(signal2FFT);
    n = pow2(nextpow2(m));  % Transform length
    
    y = fft(signal2FFT,n);          % DFT
    f = (0:n-1)*(fs/n);     % Frequency range
    
    % indices of the one-seded FFT:
    oneSideInd = 1:n/2+1;
    
    % 1-sided freq vector:
    f = f(oneSideInd);
    
    % 1-sided spectrum power vector:
    power = abs(y/n);   % 2-sided DFT
    power = power(oneSideInd);   % 1-sided DFT
    power(2:end-1) = 2*power(2:end-1); % 1-sided DFT, correct the power
    
    % 1-sided phase vector:
    phase = angle(y);
    phase = phase(oneSideInd);
    
    % find peaks in FFT in order to determine the absolute max prominence:
    [~,~,~,prominences] = findpeaks(power);
    maxProminence = max(prominences);
    
    % % % % % % find dominant peaks only:
    % freqwuency domain Process Parameters:
    min_peak_prominence = 0.05; % percent from max prominence
    min_peak_dist = 5; % the expected distance between the peaks is
                            % a multiplication of the base frequency.
    [~,dominant_freq_locs] = findpeaks(power,...
        'MinPeakProminence',min_peak_prominence*maxProminence,...
        'MinPeakDistance',min_peak_dist);

    % get fundamental harmonics:
    fundanetalHarmonics = (f(dominant_freq_locs));
    periodsFFT = 1./fundanetalHarmonics;
    
    % get the amplitude of each harmonics:
    harmonics_amp = power(dominant_freq_locs);
    
    % get the phase of each harmonics:
    harmonics_phase = phase(dominant_freq_locs);
    
    if graphicsFlag
        
        figure;
        plot(signal2FFT); grid minor;
        title('the signal on which FFT was performed:');
        
        figure;
        subplot(3,1,1)
        plot(T,X); grid minor;
        xlabel('time [sec]')
        ylabel('CPG output')
        title('CPG output over time')
        
        subplot(3,1,2)
        plot(f,power); grid minor;
        xlabel('f [Hz]');   ylabel('power');
        
        subplot(3,1,3)
        plot(f,phase(oneSideInd)*180/pi); grid minor;
        xlabel('f [Hz]');   ylabel('angle [deg]');
        
        figure;
        findpeaks(power,...
            'MinPeakProminence',min_peak_prominence*maxProminence,...
            'MinPeakDistance',min_peak_dist);
    end     

    % TODO: implement amplitude detection.
    signal_amp = [];
end

