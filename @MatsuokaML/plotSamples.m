function plotSamples( obj, results, mask, n_samples, title_str )
%PLOTSAMPLES Summary of this function goes here
%   Detailed explanation goes here

tend_temp = obj.tEnd;
tstep_temp = obj.tStep;
obj.tStep = 0.01;

if nargin<3
    n_samples = 1;
end
if nargin<4
    title_str = 'Matsuoka sample #';
end    

% Find ids of "good" samples
ids = find(mask == 1);
if isempty(ids)
    return;
end
if length(ids) < n_samples
    % Show less samples if not enough are available
    n_samples = length(ids);
end
    
if n_samples > 1
    id_samp = randsample(ids, n_samples);
else
    id_samp = ids;
end

for i = 1:n_samples
    sr = results(id_samp(i)); % Get sample by id
    samp_period = max(sr.periods);
    % Update total sim time
    if isnan(samp_period)
        obj.tEnd = tend_temp; % Use default when no period was detected
    else
        obj.tEnd = 10*max(sr.periods);
    end
    
    % Run simulation
    [~, ~, signal] = obj.runSim(sr.seq, sr.b);    
    
    % Show only the last 2 periods
    if isnan(samp_period)
        t2periods = 5;
    else
        t2periods = 2*samp_period;
    end
    t_id = find(signal.T<signal.T(end)-t2periods, 1, 'last');
    tspan = t_id:length(signal.T);
    
    figure
    hold on
    for s = 1:size(signal.signal,1)
        plot(signal.T(tspan),signal.signal(s,tspan));
    end
    % Show id+period info on title
    title([title_str,int2str(id_samp(i)),10,...
        sprintf('Period: %.3f / %.3f', sr.periods(1), sr.periods(2))]);
    % Show parameters on xlabel
    xlabel([sprintf('Params: Tr = %.3f, Ta = %.3f, b = %.3f', sr.Tr, ...
        sr.Ta, sr.b), sprintf('c = [%.3f, %.3f, %.3f, %.3f]''', ...
        sr.c(1), sr.c(2), sr.c(3),sr.c(4)), 10, ...
        sprintf('W = [%.3f, %.3f, %.3f, %.3f;', ...
        sr.W(1,1), sr.W(1,2), sr.W(1,3),sr.W(1,4)), 10, ...
        sprintf('     %.3f, %.3f, %.3f, %.3f;', ...
        sr.W(2,1), sr.W(2,2), sr.W(2,3),sr.W(2,4)), 10, ...
        sprintf('     %.3f, %.3f, %.3f, %.3f;', ...
        sr.W(3,1), sr.W(3,2), sr.W(3,3),sr.W(3,4)), 10, ...
        sprintf('     %.3f, %.3f, %.3f, %.3f]', ...
        sr.W(4,1), sr.W(4,2), sr.W(4,3),sr.W(4,4))])
    
    % Add a second figure with the full simulation result for clarity
    % (not for publishing)
    figure
    hold on
    for s = 1:size(signal.signal,1)
        plot(signal.T,signal.signal(s,:));
    end
    
%     [y, periods, signals, pos_work, neg_work] = obj.processResults(signal.X, signal.T);
end

obj.tEnd = tend_temp;
obj.tStep = tstep_temp;

