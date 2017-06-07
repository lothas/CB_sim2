function [periodsFFT,sine_coef,cos_coef,a0,fitObject] = ...
    processResults_LSQ( obj, CPG_out,T, graphicsFlag)
    % this function is for the 2N symmetric Matsuoka CPG case only
    
    % this function use 'fit' function and assume that the signal is from
    %   the shape of:
    % y = a0 + a1*cos(x*w) + b1*sin(x*w) + 
    %           a2*cos(2*x*w) + b2*sin(2*x*w) + 
    %           a3*cos(3*x*w) + b3*sin(3*x*w)
    
    % INPUTS:
    % 'T' - time vector
    % 'CPG_out' - the output of the CPG. defined as 'x2-x1'
    
    % OUTPUTS:
    % *) 'periodsFFT' - the periods of the signal (==1/harmonicsFreq)
    % *) 'sine_coef' - the coef of the sine part
    % *) 'cos_coef' - the coef of the cosine part
    % *) 'a0' - the zero freq part
    % *) 'fitObject' - the object from 'fit' function
    
    periodsFFT = NaN;
    sine_coef = NaN;
    cos_coef = NaN;
    a0 = NaN;
    fitObject = NaN;
    
    % process parameters:
    percent_of_sig2process = 0.25; % aliminate the transient part
    
    X = CPG_out;
    T = T';
    
    % Turn off findpeaks warning
    warning('off','signal:findpeaks:largeMinPeakHeight');

    dt = mean(diff(T));
    fs = 1/dt;
    
    signal_start = floor(percent_of_sig2process * length(T));
    
    signal2Process = X(1,signal_start:end);
%     % aliminate the 'zero' frequency (bias) (no need!!!!)
%     signal2Process = signal2Process - min(signal2Process);
    
    % check if the signal is non periodic
    if std(signal2Process) < 1e-3 % check standart deviation:
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
        return
    end
    % get the indices of the original signal:
    locs = signal_start+[peak_locs(1,1), ...
        peak_locs(1,2)];
    
    % prepare to perform the FFT:
    T2fit = T(1,locs(1):locs(end))';
    signal2fit = (X(1,locs(1):locs(end)))';
    
    %   curve fitting assuming furrier series:
    %   NOTE: 'fourier_i' = series of order i
    fitObject = fit(T2fit,signal2fit,'fourier3');

    
    
    % get fundamental harmonics:
    fundanetalHarmonics = [1:3]*fitObject.w;
    periodsFFT = 2*pi./fundanetalHarmonics;
    
    sine_coef = [fitObject.a1,fitObject.a2,fitObject.a3];
    cos_coef = [fitObject.b1,fitObject.b2,fitObject.b3];
    a0 = fitObject.a0;
    
    if graphicsFlag
        
        figure;
        plot(signal2fit); grid minor;
        title('the signal on which LSQ was performed:');
        
%         % % reconstruct the signal:
        w = fitObject.w;
        t = T2fit;
        y = fitObject.a0 + fitObject.a1*cos(t*w) + fitObject.b1*sin(t*w) +... 
              fitObject.a2*cos(2*t*w) + fitObject.b2*sin(2*t*w) +... 
              fitObject.a3*cos(3*t*w) + fitObject.b3*sin(3*t*w);
          
       figure;
       plot(T2fit,signal2fit,'b'); hold on
       plot(T2fit,y,'r'); grid minor;
       xlabel('time [sec]');    ylabel('signal');
       legend('original signal','reconstruted signal using LSQ');
       title('comparison of fitting');
       
    end     

    % TODO: implement amplitude detection.
    signal_amp = [];
end

