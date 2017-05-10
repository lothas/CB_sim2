function [percent_osc_new,conv_in_range,accuracy] = ...
    NN_GA_perf(net,sampl,results_old,fileName,seqOrder,caseNum,perORfreq)
% this function change a sript to automate the NN perf checking using GA
% FOR the 4N Matsuoka general CP only!

% Inputs:
% *) ' net' - the NN object
% *) 'sampl' - NN inputs
% *) 'results_old' - the data of the genes before change
% *) 'fileName' - a file name to the sim results (small file please!)
% *) 'seqOrder' - the gene seqence
% *) 'caseNum' - which case are we checking
% *) 'perORfreq' - are we using the period or the frequency

% Outputs:
% *) 'percent_osc_new' - home many new CPG's are oscillating 
% *) 'conv_in_range' - home many new CPG's are oscillating in the desired
%                       range
% *) 'accuracy' - look at accuracy from the paper


%% % % % 0 stage - initiate MML:
% Initialize machine learning object for Matsuoka analysis
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01; % 0.01
MML.tEnd = 30;
nPlotSamples = 0; % 10;
% Turn off findpeaks warning
warning('off','signal:findpeaks:largeMinPeakHeight');

%% % % % 1st stage: pass the input to the NN:

theta_S1_new = net(sampl);

switch caseNum
    case {1,2,3,4}
        theta_S1_new_MAX = 0.25; % NOTE: update accordinly to "'MatsuokaGenome.mat'" file
    case {5,6,7,8}
        theta_S1_new_MAX = 10;
    case {9,10}
        theta_S1_new_MAX = [0.25;10];
    otherwise
        error('illigal caseNUM');
end
% can't select NN output which is out of bound:
good_ids = (theta_S1_new > 0.02) & (theta_S1_new < theta_S1_new_MAX);
if size(theta_S1_new,1) == 2
    good_ids = good_ids(1,:) & good_ids(2,:);
end
results_old = results_old(good_ids);
theta_S1_new = theta_S1_new(:,good_ids);

%% % % % 2nd stage: - make new CPGs:
[results_new] = ...
    change_CPG(MML,caseNum,seqOrder,theta_S1_new,results_old);

%% % % % 3rd stage: results:
periods_new = horzcat(results_new(:).periods);
periods_new = periods_new(1,:);
% osc_ids_new = find(~isnan(periods_new));

percent_osc_new = (sum(~isnan(periods_new))/length(periods_new));
disp(['the percentage of osc period after the change is: ',...
    num2str(100*percent_osc_new),'%']);

conv_in_range_temp = (periods_new > MML.perLimOut(1,1)) & ...
    (periods_new < MML.perLimOut(1,2));
conv_in_range = sum(conv_in_range_temp) / length(periods_new);
disp(['the percentage of CPGs which converge in period range is: ',...
    num2str(100*conv_in_range),'%']);

periods_old = horzcat(results_old(:).periods);
periods_old = periods_old(1,:);
delta_vec = ((1./periods_old) .* ( periods_new - periods_old ));
delta = mean(delta_vec,'omitnan');
accuracy = 1/(1+delta);

disp(['the Accuracy is: ',...
    num2str(accuracy)]);
end