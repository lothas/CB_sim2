% picks = randsample(24000, 50);
% picks = [24893, 10746, 2447, 6292, 3632, 9959];
picks = [1408, 15091, 14144, 8277, 27433, 27893, 10406, 15161, 24716, 26680];

for i = 1:length(picks)
    MML.tEnd = 5;
    MML.tStep = 0.02;
    sr = results(picks(i)); % Get sample by id
    [out, ~, signal] = MML.runSim(sr.seq, sr.b);
    disp(out.periods')
    
    figure
    hold on
    plot(signal.T,signal.signal);
    % Show id+period info on title
    title(['Sample #',int2str(picks(i)),10,...
        sprintf('Period: %.3f / %.3f', out.periods(1), out.periods(2))]);
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
    
    if 0
        testing1 = 1;
        [y, periods, signals, pos_work, neg_work] = MML.processResults(signal.X, signal.T);
        [simFreq, amp] = MML.processResultsFFT(signal.X, signal.T, 1);
        testing1 = 0;
    end
end