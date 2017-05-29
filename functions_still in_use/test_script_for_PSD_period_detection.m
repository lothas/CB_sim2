
close all; clear all; clc

%% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 50; % 15
MML.nNeurons = 2;
%%
filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT.mat';
load(filename,'results');

%%
% picks = randsample(length(results), 3);
% picks = [24893, 10746, 2447, 6292, 3632, 9959];
picks = 2;

for i = 1:length(picks)
    sr = results(picks(i)); % Get sample by id
    [out, ~, signal] = MML.runSim(sr.seq);%, sr.b);
%     seq = sr.seq;
%     [out, ~, signal] = MML.runSim(seq);
    disp(['the output period from autoCorr is ',...
        num2str(out.periods'),' [seq]']);
    
    % For 2N symm CPG case:
    if 1
        figure
        hold on
        plot(signal.T,signal.signal);
        % Show id+period info on title
        title(['Sample #',int2str(picks(i)),10,...
            sprintf('Period: %.3f ', out.periods(1))]);
        % Show parameters on xlabel
        % % % % % % ASUMMING: T = 5*tau !!; % % % % % %
        xlabel([sprintf('Params: Tr = %.3f, Ta = %.3f, b = %.3f', sr.tau, ...
            5*sr.tau, sr.b), ...
            sprintf('c = %.3f', sr.c(1)), ...
            sprintf('a = %.3f;', sr.a)])
    end
    
    % For 4N CPG case:
    if 0
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
    end
    
end

%% PSD of the signal:
MML.processResultsFFT_better(signal.X, signal.T, 1);