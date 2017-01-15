

clc; close all; clear all;

%% Load varibles:
% load('MatsRandomRes_test.mat')

load('MatsRandomRes_16_12_2016.mat','results','periods')
results1 = results; periods1 = periods;
clear results periods
load('MatsRandomRes_18_12_2016.mat','results','periods')
results2 = results; periods2 = periods;
clear results periods
load('MatsRandomRes_19_12_2016.mat','results','periods')
results3 = results; periods3 = periods;
clear results periods

results = horzcat(results1,results2,results3);
periods = horzcat(periods1,periods2,periods3);

clear results1 results2 results3
clear periods1 periods2 periods3
%% 
HiddenN = 10;
NumOfRepeats = 20;
netMseTrain = zeros(NumOfRepeats,8);
netMseValidation = zeros(NumOfRepeats,8);
netMseTest = zeros(NumOfRepeats,8);

ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);
clear ids_period ids_error nSims

targetCells = {'freq'};

for i=1:8
    switch i
        case 1
            parametersCells = {'tau'};
            W_flag = false;
        case 2
            parametersCells = {'tau','b'};
            W_flag = false;
        case 3
            parametersCells = {'tau','b',...
                                'w_{12}','w_{13}','w_{14}',...
                                'w_{21}','w_{23}','w_{24}',...
                                'w_{31}','w_{32}','w_{34}',...
                                'w_{41}','w_{42}','w_{43}'};
            W_flag = false;
        case 4
            parametersCells = {'tau','b',...
                                'w_{12}','w_{13}','w_{14}',...
                                'w_{21}','w_{23}','w_{24}',...
                                'w_{31}','w_{32}','w_{34}',...
                                'w_{41}','w_{42}','w_{43}'};
            W_flag = true;
        case 5
            parametersCells = {'tau','b','prodW','sumW'};
            W_flag = false;
        case 6
            parametersCells = {'tau','b','prodW','sumW'};
            W_flag = true;
        case 7
            parametersCells = {'tau','b','W123','W124','W132','W134',...
                'W142','W143','W234','W243'};
            W_flag = false;
        case 8
            parametersCells = {'tau','b','W123','W124','W132','W134',...
                'W142','W143','W234','W243'};
            W_flag = true;
        otherwise
            error('i not in range');
    end
    
    [ sampl,targ ] = prepareData( results,periods,ids,parametersCells,targetCells,W_flag );
    
    normParams = zeros(size(sampl, 1), 2); % Normalize samples
    for j = 1:size(sampl, 1)
        feat = sampl(j, :);
        normParams(j, :) = [mean(feat), std(feat)];
        sampl(j, :) = (feat - normParams(j, 1))/normParams(j, 2);
    end
    clear feat normParams
    
    for j=1:NumOfRepeats
        net = feedforwardnet(HiddenN);
        [net, tr] = train(net, sampl, targ);
        netMseTrain(j,i) = tr.best_perf;
        netMseValidation(j,i) = tr.best_vperf;
        netMseTest(j,i) = tr.best_tperf;
        clear net tr
    end
end

% save('NNdata.mat','netMseTrain','netMseValidation','netMseTest','HiddenN')
meanMseTrain = mean(netMseTrain,1);
meanMseValidation = mean(netMseValidation,1);
meanMseTest = mean(netMseTest,1);
means = [meanMseTrain;meanMseValidation;meanMseTest];

stdMseTrain = std(netMseTrain,0,1);
stdMseValidation = std(netMseValidation,0,1);
stdMseTest = std(netMseTest,0,1);
stdevs = [stdMseTrain;stdMseValidation;stdMseTest];


Names={' ';'tau                       ';'tau,b                             ';
       'tau,b,w_{ij}              ';'tau,b,w_{ij} hat                  ';
       'tau,b,Pi w_{ij},sum w_{ij}';'tau,b,Pi w_{ij} hat,sum w_{ij} hat';
       'tau,b,w_{ijk}             ';'tau,b,w_{ijk} hat                 '};

numgroups = size(means, 2); 
numbars = size(means, 1); 
groupwidth = min(0.8, numbars/(numbars+1.5));

figure; hold on
h=bar(means');
set(gca,'XTickLabel',Names);
set(h,'BarWidth',1);
for k=1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*k-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, means(k,:), stdevs(k,:), 'k', 'linestyle', 'none');
end
legend('train','valid','test');
ylabel('MSE');
title('MSE over different NN (target=frequency)');
hold off

% figure; hold on
% bar(1:8,meanMseTrain);
% errorbar(1:8,meanMseTrain,stdMseTrain,'.');
% set(gca,'XTickLabel',Names);
% title('MSE over different NN (training)(target=frequency)');
% ylabel('MSE');
% 
% figure; hold on
% bar(1:8,meanMseValidation);
% errorbar(1:8,meanMseValidation,stdMseValidation,'.');
% set(gca,'XTickLabel',Names);
% title('MSE over different NN (valid)(target=frequency)');
% ylabel('MSE');
% 
% figure; hold on
% bar(1:8,meanMseTest);
% errorbar(1:8,meanMseTest,stdMseTest,'.');
% set(gca,'XTickLabel',Names);
% title('MSE over different NN (test)(target=frequency)');
% ylabel('MSE');

%% Statistical test
x = netMseTest(:,3);
y = netMseTest(:,7);
[h_ttest,p_ttest,ci_ttest,stats_ttest] = ttest(x,y);
[p_stest,h_stest,stats_stest] = signtest(x,y);


