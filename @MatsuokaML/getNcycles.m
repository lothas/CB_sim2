function [sigT, signal] = getNcycles(obj, T, in_signal, period, N)
%GETNCYCLES Summary of this function goes here
%   Detailed explanation goes here
    if any(isnan(period))
        sigT = T;
        signal = in_signal;
        return;
    end
        
    last_id = [];
    pos = length(T)-1;
    
    if min(in_signal(1,floor(pos/2):end))>0 || ...
            max(in_signal(1,floor(pos/2):end))<0
        % Move signal to x axis
        in_signal = bsxfun(@minus, in_signal, mean(in_signal,2));
    end
    
    % Find last "crossing" pointpos = length(T)-1;
    while pos>0
        if xor(in_signal(1,pos)>0,in_signal(1,pos+1)>0)
            % zero crossing
            if (in_signal(1,pos)<in_signal(1,pos+1))
                % crossing from negative to positive
                last_id = pos;
                break;
            end
        end
        pos = pos-1;
    end
    
    if isempty(last_id)
        % Still couldn't find a crossing point (shouldn't happen)
        sigT = T;
        signal = in_signal;
        return;
    end
    
    dt = mean(diff(T));
    first_id = max(1, floor(last_id - N*max(period)/dt));
    signal = in_signal(:,first_id:last_id);
    
    sigT = (1:size(signal,2))*dt;
end

