clear all; close all; clc

%% Load Data:
results1 = load('MatsRandomRes_4Neurons_with_LSQ_1.mat','results');
results2 = load('MatsRandomRes_4Neurons_with_LSQ_2.mat','results');

results = horzcat(results1.results,results2.results);
clear results1 results2

%% get periods:
periods_AC = horzcat(results(:).periods);
periods_LS = (vertcat(results(:).periods_LSQ))';

periods_AC_hip = periods_AC(1,:);
periods_AC_ankle = periods_AC(2,:);

periods_LS_hip = periods_LS(1,:);
periods_LS_ankle = periods_LS(2,:);

periods_AC_mean = mean(periods_AC,1);
periods_LS_mean = mean(periods_LS,1);
% TODO: do something is one of the periods os NaN
%% phase 1- confusion matrix, osc/non-osc:
% for ankle torque:
ids_AC_ankle = ~isnan(periods_AC_ankle);
ids_LS_ankle = ~isnan(periods_LS_ankle);

figure;
plotconfusion(ids_AC_ankle,ids_LS_ankle,'Ankle torque')
xlabel('Auto-Corr')
set(gca,'xticklabel',{'n-osc' 'osc' ''})
ylabel('LS')
set(gca,'yticklabel',{'n-osc' 'osc' ''})

% for ankle torque:
ids_AC_hip = ~isnan(periods_AC_hip);
ids_LS_hip = ~isnan(periods_LS_hip);

figure;
plotconfusion(ids_AC_hip,ids_LS_hip,'Hip torque')
xlabel('Auto-Corr')
set(gca,'xticklabel',{'n-osc' 'osc' ''})
ylabel('LS')
set(gca,'yticklabel',{'n-osc' 'osc' ''})

%% correlation coefficient between AC and LS
best_ids_ankle = ids_AC_ankle & ids_LS_ankle;
best_ids_hip = ids_AC_hip & ids_LS_hip;

R_ankle = corrcoef(periods_AC_ankle(1,best_ids_ankle),...
    periods_LS_ankle(1,best_ids_ankle))

R_hip = corrcoef(periods_AC_hip(1,best_ids_hip),...
    periods_LS_hip(1,best_ids_hip))
%% repeat with mean periods
ids_AC_mean = ~isnan(periods_AC_mean);
ids_LS_mean = ~isnan(periods_LS_mean);

figure;
plotconfusion(ids_AC_mean,ids_LS_mean,'Mean torque')
xlabel('Auto-Corr')
set(gca,'xticklabel',{'n-osc' 'osc' ''})
ylabel('LS')
set(gca,'yticklabel',{'n-osc' 'osc' ''})

%% check conditions from section 3.5 in paper:
y_i_osc = (vertcat(results(:).neuronOsc))';

y_12_osc = y_i_osc(1,:) | y_i_osc(2,:);
y_34_osc = y_i_osc(3,:) | y_i_osc(4,:);

y_12_and_34_osc = y_12_osc & y_34_osc;

periods_AC = horzcat(results(:).periods);
osc_net = ~isnan(periods_AC(1,:)) & ~isnan(periods_AC(2,:));
%% run all simulations again and check the conditions in sec 3.5 in the paper:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 30; % 15
MML.nNeurons = 4;

% N = length(results);
N=100;

disp('start with the sim:');
parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
% for i=1:N
    disp(['at sim #',num2str(i)]);
    
    period_peaks = nan(1,2);
    sq_err = nan(1,2);
    firstCheck = false(1,2);
    
    [out, sim, signals] = MML.runSim(results(i).seq);
        % Prepare output:
    % Parameters
    results_more(i).seq = results(i).seq;

    % Results- caculate perdiods using different methods:
    results_more(i).periods = out.periods;
    
    [periods_LSQ1,~,~,~,~] = ...
        MML.processResults_LSQ(signals.signal(1,:),signals.T,0);
    [periods_LSQ2,~,~,~,~] = ...
        MML.processResults_LSQ(signals.signal(2,:),signals.T,0);
    results_more(i).periods_LSQ = [periods_LSQ1(1,1),periods_LSQ2(1,1)];
    
    if ~isnan(out.periods(1,1)) && ~isnan(out.periods(2,1))
        % I) Verify that the squared errors between T_out and the distance
        %   between the first and second peaks of each output signal 
        %   are below 1E-4.
        sig2proc = 0.3; % start calculating period after 30% of the signal
                        % to skip transient
        sig_idx = floor(sig2proc*size(signals.signal,2)):size(signals.signal,2);

        for j=1:2
            signal = signals.signal(j,:);
            norm_signal = (signal(sig_idx)-min(signal(sig_idx))) / ...
                        (max(signal(sig_idx)) - min(signal(sig_idx)));
            [~,locs]=findpeaks(norm_signal, 'MinPeakheight',0.5, ...
                'MinPeakProminence', 0.05);

            if length(locs) > 3
                period_peaks(1,j) = signals.T(locs(2)) - signals.T(locs(1));
                sq_err(1,j) = (period_peaks(1,j) - out.periods(j,1))^2;
                firstCheck(1,j) = (sq_err(1,j) < 1e-4);
            end
        end
    
    end
    results_more(i).period_peaks = period_peaks;
    results_more(i).sq_err = sq_err;
    results_more(i).firstCheck = firstCheck;
    
    results_more(i).neuronActive = results(i).neuronActive;
    results_more(i).neuronOsc = results(i).neuronOsc;
    
    
    % II) Verify that y_1 or y_2 and y_3 or y_4 were oscillatory, detected
    %   as sigma(y_i)>1E-4 , if a period was detected and not oscillatory
    %   otherwise
    
end 
disp('sim end...');

save('secnd_varification.mat','results');