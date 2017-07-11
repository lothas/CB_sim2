
clc; close all; clear all;

%%
filename = 'MatsRandomRes_2Neurons_symm_test_for_shapeLearn_changing_C.mat';

load(filename,'results');

periods_AC = vertcat(results(:).periods);
periods_LSQ = zeros(size(periods_AC));

for i=1:300000
    temp1 = results(i).periods_LSQ;
    periods_LSQ(i,1) = temp1(1,1);
end

%% plot the Histograms of the two to check matches:
figure;
histogram(periods_AC,100,'FaceAlpha',0.5); hold on;
histogram(periods_LSQ,100,'FaceAlpha',0.2);
legend('autoCorelation','NL-LSQ');

%%
AC_ids = ~isnan(periods_AC);
LSQ_ids = ~isnan(periods_LSQ);
disp(['osc from autoCorr is:   ',num2str(sum(AC_ids))]);
disp(['osc from LSQ is:   ',num2str(sum(LSQ_ids))]);

ids = AC_ids & LSQ_ids;
disp(['osc from both is:   ',num2str(sum(ids))]);
%% plot Histogram of the error between the two methods:
Per_err = periods_AC(:,1) - periods_LSQ(:,1);
figure;
edges = [-10 -0.5:0.005:0.5 10];
histogram(Per_err,edges)