function [fundanetalHarmonics, amp] = processResults_PSD( obj, X, T, graphicsFlag)
    %checking the period using limit cycle 
    
    % INPUTS:
    % T - time vector
    % X - the state space
    % samplingTime - the sampling time of ode45
    % ssTime - the time which we estimate arrival to steady state (limit
    %          cycle).
    
    X = X';
    T = T';
    % Turn off findpeaks warning
    warning('off','signal:findpeaks:largeMinPeakHeight');

    % Process Parameters:
    min_peak_prominence = 0.7; % percent from max prominence, freq domain
    min_peak_hight = 0.85; % percent from max hight, time domain
    samplingTime = mean(diff(T));
    Fs = 1/samplingTime;
%     NFFT = 2^nextpow2(length(T));

    % Build signals
    if obj.nNeurons > 2
        y = max(X(1:obj.nNeurons,:), 0)';
        signals = obj.Sim.Con.OutM*y;
    elseif obj.nNeurons == 2
        % NOTE: if nNuerons=2 than there is always a zero vector because of
        %   M_out (row of zeros).
        y = max(X(1:2,:), 0);
        signals = [1,-1]*y;
    else
        error('invalid number of neurons')
    end
    
    nSignals = size(signals, 1);
    peak_locs = zeros(nSignals,2);
    signal_start = floor(0.25*size(signals,2));
    longest_signal = 1; signal_length = 1;
    for s = 1:nSignals
        signal2Process = signals(s,signal_start:end);
        signal2Process = signal2Process - min(signal2Process);
        
        % check if the signal is non oscillatory. method #1
        if std(signal2Process) < 1e-3
            fundanetalHarmonics = NaN;
            % TODO: implement amplitude detection
            amp = []; % mean(signals'); %#ok<UDIM>
            return
        end
        
        [~,locs]=findpeaks(signal2Process, ...
            'MinPeakheight',...
            min_peak_hight*max(signal2Process));

        if length(locs) > 1
            peak_locs(s,:) = [locs(1), locs(end)];
        else
            fundanetalHarmonics = NaN;
            amp = mean(signals'); %#ok<UDIM>
            return
        end
        
        % ????
        if (peak_locs(s,2)-peak_locs(s,1) > signal_length)
            signal_length = peak_locs(s,2)-peak_locs(s,1);
            longest_signal = s;
        end
    end
     
    locs = signal_start+[peak_locs(longest_signal,1), ...
        peak_locs(longest_signal,2)];
    
    signal = signals(:,locs(1):locs(end));
    T2graph = T(1,locs(1):locs(end));
    
    % PSD parameters:
    window = hamming(512);
    
    % PSD:
    try
        [Pxx,Fxx] = pwelch(signal,window,[],[],Fs);
    catch
        % dump oscillating signals get falsly flaged as oscillating.
        %   which leads to "signals" with really small num of samples. 
        warning('was unable to compute "pwelch"')
        fundanetalHarmonics = NaN;
        % TODO: implement amplitude detection
        amp = [];
        return
    end
%         [Pxx,Fxx] = pwelch(x,window,noverlap,f,fs)
    Pxx_dB = 10*log10(Pxx);
    
    % find peaks in PSD:
    [~,~,~,prominences] = findpeaks(Pxx_dB);
    maxProminence = max(prominences);

    % find dominant peaks only:
    [~,locs] = findpeaks(Pxx_dB,'MinPeakProminence',...
        min_peak_prominence*maxProminence);
    
    fundanetalHarmonics = (Fxx(locs))';
    
    if graphicsFlag
        figure;
        subplot(2,1,1);
        plot(Fxx,Pxx_dB); grid minor;
        xlabel('Hz');
        ylabel('dB')

        subplot(2,1,2);
        plot(T2graph,signal); grid minor;
        xlabel('Time [sec]');
        ylabel('signal')

        figure; % display peaks in a sorted way:
        [psor,lsor] = findpeaks(Pxx_dB,'SortStr','descend');
        findpeaks(Pxx_dB);
        title('sorted peaks big to small');
        text(lsor+.02,psor,num2str((1:numel(psor))'))
        clear psor lsor
        
        figure;
        findpeaks(Pxx_dB,'MinPeakProminence',...
            min_peak_prominence*maxProminence)
        title('peaks with most prominence:');
        xlabel('locations');    ylabel('Magnitude [dB]');
        
%         % reconstruct the signal:
%         harmonics_amp = Pxx(locs);
% %         harmonics_amp = sqrt(harmonics_amp*Fs/2)
%         signal_from_PSD = harmonics_amp(1,1) *...
%             sin(2*pi*fundanetalHarmonics(1,1)*T);
%         if length(fundanetalHarmonics) > 1
%             for j=2:length(fundanetalHarmonics)
%                 signal_from_PSD = signal_from_PSD + (harmonics_amp(j,1) *...
%                 sin(2*pi*fundanetalHarmonics(1,j)*T));
%             end
%         end
%         figure;
%         plot(T,signal_from_PSD); hold on
%         plot(T,signals);
%         xlabel('time'); ylabel('signal');
%         legend('FFT reconstructed signal','original signal');
    end     

    % TODO: implement amplitude detection.
    amp = [];
end

