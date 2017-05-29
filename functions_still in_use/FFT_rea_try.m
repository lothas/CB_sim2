
close all; clear all; clc

%% define signal:
N = 2000;
tmin = 0;
tmax = 30;
T = linspace(tmin,tmax,N);

sampling_rate = N/(tmax-tmin);
% % simple sine wave:
% x = 3*sin(2*pi*5*T);
% % sine wave with random noise:
% x = 3*sin(2*pi*5*T) + 5*sin(2*pi*16*T)+rand(1,length(T));
% % simple aquare wave:
x = 0.5 + 0.5 * square(2*pi*5*T,50);
%% FFT
% INPUTS:
    % T - time vector
    % X - the state space
    % samplingTime - the sampling time of ode45
    % ssTime - the time which we estimate arrival to steady state (limit
    %          cycle).


    
signals = x;    
samplingTime = mean(diff(T));
Fs = 1/samplingTime;

signal_start = floor(0.25*size(x,2));
longest_signal = 1; signal_length = 1;
nSignals = size(signals, 1);
peak_locs = zeros(nSignals,2);

for s = 1:nSignals
    signal2Process = signals(s,signal_start:end);
    % bring the signal above the x axis:
    signal2Process = signal2Process - min(signal2Process);

    % don't bother with non oscillatory signals:
    if std(signal2Process) < 1e-3
        simFreq = NaN;
        amp = mean(signals'); %#ok<UDIM>
        return
    end

    [~,locs]=findpeaks(signal2Process, ...
        'MinPeakheight',0.75*max(signal2Process));

    if length(locs) > 1
        peak_locs(s,:) = [locs(1), locs(end)];
    else
        simFreq = NaN;
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

% PSD parameters:
window = hamming(1024);

% PSD:
[Pxx,Fxx] = pwelch(signal,window,[],[],Fs);
%     [Pxx,Fxx] = pwelch(x,window,noverlap,f,fs)
Pxx_dB = 10*log10(Pxx);

if 1 %graphicsFlag
    figure;
    subplot(2,1,1);
    plot(Fxx,Pxx_dB); grid minor;
    xlabel('Hz');
    ylabel('dB')

    subplot(2,1,2);
    plot(T(1,locs(1):locs(end)),signal); grid minor;
    xlabel('Time [sec]');
    ylabel('signal')
    
    figure; % display peaks in a sorted way:
    [psor,lsor] = findpeaks(Pxx_dB,'SortStr','descend');
    findpeaks(Pxx_dB);
    text(lsor+.02,psor,num2str((1:numel(psor))'))
    clear psor lsor
end     
 
% find peaks in PSD:
[~,~,~,prominences] = findpeaks(Pxx_dB);
maxProminence = max(prominences);

[~,locs] = findpeaks(Pxx_dB,'MinPeakProminence',maxProminence/2);
fundanetalHarmonics = Fxx(locs);

figure;
findpeaks(Pxx_dB,'MinPeakProminence',maxProminence/2)
title('peaks with most prominence:');
xlabel('locations');    ylabel('Magnitude [dB]');