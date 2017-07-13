function [percent_osc_new,conv_in_range,accuracy] = ...
    NN_GA_perf(net,sampl,results_old,seqOrder,caseNum)
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
numOfCPGs = length(theta_S1_new);
% NOTE: update accordinly to "'MatsuokaGenome.mat'" file!!
tau_min = 0.02;
tau_Max = 0.25;
b_min = 0.2;
b_max = 10;

% can't select NN output which is out of bound:
switch caseNum
    case {1,2,3,4}
        good_ids = (theta_S1_new > tau_min) & ...
                (theta_S1_new < tau_Max); 
    case {5,6,7,8}
        good_ids = (theta_S1_new > b_min) & ...
                (theta_S1_new < b_max);
    case {9,10}
        good_ids = (theta_S1_new(1,:) > tau_min) & ...
            (theta_S1_new(1,:) < tau_Max) &...
            (theta_S1_new(2,:) > b_min) & ...
            (theta_S1_new(2,:) < b_max);
    otherwise
        error('illigal caseNUM');
end

results_old = results_old(good_ids);
theta_S1_new = theta_S1_new(:,good_ids);

disp(['the number of CPGs with parameters in bounds is ',...
    num2str(sum(good_ids)),' out of ',num2str(numOfCPGs)]);
%% % % % 2nd stage: - make new CPGs:
[results_new] = ...
    change_CPG(MML,caseNum,seqOrder,theta_S1_new,results_old);

%% extract periods before NN (desired periods) and after (actual periods)
periods_old = horzcat(results_old(:).periods);
periods_old = periods_old(1,:);
periods_new = horzcat(results_new(:).periods);
periods_new = periods_new(1,:);
%% Plot periods histograms of before and after the change
% figure;
% histogram(periods_old,100,'Normalization','pdf'); hold on;
% histogram(periods_new,100,'Normalization','pdf');
% legend('before NN','after NN');
% xlabel('periods [sec]');
%% % % % 3rd stage: results:

percent_osc_new = ( sum(~isnan(periods_new)) / length(periods_new) );
disp(['the percentage of osc period after the change is: ',...
    num2str(100*percent_osc_new),'%']);

conv_in_range_temp = (periods_new > MML.perLimOut(1,1)) & ...
    (periods_new < MML.perLimOut(1,2));
conv_in_range = sum(conv_in_range_temp) / length(periods_new);
disp(['the percentage of CPGs which converge in period range is: ',...
    num2str(100*conv_in_range),'%']);

periods_des = (MML.perLimOut(1,2)+MML.perLimOut(1,1))/2;
delta_vec = ((1./periods_des) .* ( periods_new - periods_des ));
% delta_vec = (( periods_new - periods_old )./periods_old);
delta = mean(delta_vec,'omitnan');
accuracy = 1/(1+delta);

disp(['the Accuracy is: ',...
    num2str(accuracy)]);
end