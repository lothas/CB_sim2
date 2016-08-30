function [y, periods, signals, pos_work, neg_work] = ...
        processResults(obj, X, T)
%PROCESSRESULTS Calculates the output signals, CPG period, positive and
%negative work using the simulation results
    y = max(0, X(:, 1:2:end));
    % Initialize output variables
    periods = zeros(obj.nNeurons/2,1);
    signals = zeros(length(y), obj.nNeurons/2);
    pos_work = zeros(obj.nNeurons/2,1);
    neg_work = zeros(obj.nNeurons/2,1);

    for i = 1:obj.nNeurons/2
        signal = y(:,2*i-1)-y(:,2*i);

        % Calculate signal period
        sig2proc = 0.3; % start calculating period after 30% of the signal
        % to skip transient
        sig_idx = floor(sig2proc*length(signal)):length(signal);
        ac = xcorr(signal(sig_idx),signal(sig_idx),'coeff');
        [~,locs]=findpeaks(ac, 'MinPeakheight',0.3);
        periods(i) = mean(diff(locs))*mean(diff(T));

        % Calculate positive and negative work
        signals(:,i) = signal;
        pos_signal = max(signal,0);
        neg_signal = min(signal,0);
        pos_work(i) = trapz(T,pos_signal);
        neg_work(i) = trapz(T,neg_signal);
    end
end

