
close all; clear all; clc

% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 50; % 15
MML.nNeurons = 2;
%
% filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT.mat';
filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT.mat';
load(filename,'results');
%% Run individual simulations:
%
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

x = signal.signal(1,:);
t = signal.T;

%% LSQ
close all; clc;

[periodsFFT,sine_coef,cos_coef,a0,fitObject] = ...
    MML.processResults_LSQ(x,t,1)

close all
N = 10000;
t=linspace(0,3,N);
x = sin(100*t) + 5*rand(1,length(t));
MML.processResults_FFT(x,t,1)
%% FFT
close all; clc;

[periodsFFT,f_amp,f_phase, ~] = ...
    MML.processResults_FFT(x,t,1)

f_har = 1./periodsFFT;

x_rec = 0;
for i=1:length(f_har)
    x_rec = x_rec + ...
        (f_amp(1,i) .* sin(2*pi*f_har(1,i).*t + f_phase(1,i)));
end

figure;
plot(t,x,'b'); hold on;
plot(t,x_rec,'r'); grid minor;
xlabel('time [sec]');
ylabel('signal');
title(' comparison between the actual signal and the reconstruction with FFT');
legend('original signal','from FFT');


%% checking known signals
close all; clc;

N = 5100;
t = linspace(0,20,N);

% % beating with noise:
x = 0.5*sin(2*pi*10*t) + 0.3*sin(2*pi*11*t) + 0.2*rand(1,length(t));

% % % simple sine:
% x = 0.5*sin(2*pi*10*t);

% % % squar wave:
% x = 0.5 + 0.5 * square(2*pi*t,50);

