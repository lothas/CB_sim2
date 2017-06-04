function [simFreq, amp] = processResultsFFT( obj, X, T, graphicsFlag)
    %checking the period using limit cycle 
    
    % INPUTS:
    % T - time vector
    % X - the state space
    % samplingTime - the sampling time of ode45
    % ssTime - the time which we estimate arrival to steady state (limit
    %          cycle).
    
    % Turn off findpeaks warning
    warning('off','signal:findpeaks:largeMinPeakHeight');

    samplingTime = mean(diff(T));
    
    % FFT
%     smallestFrequency = 1/(t(end)-ssTime);

    % Build signals
    y = max(X(:,1:obj.nNeurons), 0)';
    signals = obj.Sim.Con.OutM*y;
    
    nSignals = size(signals, 1);
    peak_locs = zeros(nSignals,2);
    signal_start = floor(0.25*size(signals,2));
    longest_signal = 1; signal_length = 1;
    for s = 1:nSignals
        signal2Process = signals(s,signal_start:end);
        signal2Process = signal2Process - min(signal2Process);
        
        if std(signal2Process) < 1e-3
            simFreq = NaN;
            amp = mean(signals'); %#ok<UDIM>
            return
        end
        
        [~,locs]=findpeaks(signal2Process, ...
            'MinPeakheight',0.75*max(signal2Process));
    %     [~,locs]=findpeaks(signal2Process, 'Threshold',abs(max(signalTemp)/10));

        if length(locs) > 1
            if length(locs) > 9
                peak_locs(s,:) = [locs(end-9), locs(end)];
            else
                peak_locs(s,:) = [locs(1), locs(end)];
            end
        else
            simFreq = NaN;
            amp = mean(signals'); %#ok<UDIM>
            return
        end
        
        if (peak_locs(s,2)-peak_locs(s,1) > signal_length)
            signal_length = peak_locs(s,2)-peak_locs(s,1);
            longest_signal = s;
        end
    end
     
    locs = signal_start+[peak_locs(longest_signal,1), ...
        peak_locs(longest_signal,2)];
    
    signal = signals(:,locs(1):locs(end));

%     % Reality check
%     T = T(locs(1):locs(end));
%     for s = 1:nSignals
%         signal2Process = signal(s,:);
%         signal2Process = signal2Process - min(signal2Process);
%         [~,locs]=findpeaks(signal2Process, ...
%             'MinPeakheight',0.8*max(signal2Process));
%         disp(mean(diff(T(locs))))
%     end
    
    Fs = 1/samplingTime;                    % Sampling frequency
    L = length(signal);                     % Length of signal
    n = 2^nextpow2(L);
    Y = fft(signal',n);

    P2 = abs(Y/n);
    P1 = P2(1:n/2+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);

    FreqVec = 0:(Fs/n):(Fs/2-Fs/n);
    oneSidedFFT = P1(1:n/2,:);
%     plot(FreqVec,oneSidedFFT)
%     plot(1./FreqVec,oneSidedFFT)
%         oneSidedFFT_smooth = smooth(oneSidedFFT,10);

    HowManyHarmonics = 5; %how many harmonics to look for
    nFirstHarmonics = zeros(1,HowManyHarmonics);
    k=1;
    neighbors = floor(0.01*(n/2)); %how many neighboring points to check
    for j=(neighbors+1):(length(oneSidedFFT)-neighbors)
        if k==(HowManyHarmonics+1)
            break; % stop the for loop when we find the necessary harmonics
        end

        % checking if its a truly local maxima = actual harmonic
        Vec = oneSidedFFT((j-neighbors):(j+neighbors),:);
        if any(max(Vec) == Vec(neighbors+1,:))
            nFirstHarmonics(k) = j;
            k = k + 1;
        end
    end

    dominantHarmonicIndex = nFirstHarmonics(1);
    % TODO: I might have problem if the 1st harmonic is really close to
    % the zero harmonic

    if graphicsFlag
        figure;
        plot(FreqVec,oneSidedFFT);
        title('in the Frequency Domain');
        xlabel('Frequency [Hz]'); ylabel('power');

%             figure;
%             findpeaks(oneSidedFFT, 'MinPeakheight',0.3*max(oneSidedFFT))
%             findpeaks(oneSidedFFT_smooth);

    end

    if dominantHarmonicIndex
        simFreq = FreqVec(dominantHarmonicIndex);
    else
        simFreq = NaN;
    end
            
    % check amplitude:
%     amp = (max(signal2Process-mean(signalTemp))-min(signal2Process-mean(signal2Process)))/2;
    amp = max(max(signal,[],2)-mean(signal,2), ...
                mean(signal,2)-min(signal,[],2));
    % TODO: correct the amp calc.
    
    
end

