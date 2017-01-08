

clc; close all; clear all;

%% Load varibles:
load('MatsRandomRes_16_12_2016.mat','nSims','results','periods')


%% 
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);

% parametersCells = {'tau','b',...
%     'w_{12}','w_{13}','w_{14}',...
%     'w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}',...
%     'w_{41}','w_{42}','w_{43}'};
% parametersCells = {'tau','b',...
%     'prodW','sumW'};
% parametersCells = {'tau','b'};
% parametersCells = {'tau'};
parametersCells = {'tau','b',...
    'W123','W124','W132','W134','W142','W143','W234','W243'};

targetCells = {'freq'};
[ sampl,targ ] = prepareData( results,periods,ids,parametersCells,targetCells,0 );

figure;
histogram(max(horzcat(results(ids).perError2)',[],2),1000);
title('Hist of period error'); xlabel('period error');

clear ids_period ids_error nSims
clear results periods

%% Normalize samples
normParams = zeros(size(sampl, 1), 2);
for i = 1:size(sampl, 1)
    feat = sampl(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    sampl(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    
end

% targnotnorm = targ;
% normParams_targ = zeros(size(sampl, 1), 2);
% for i = 1:size(targ, 1)
%     feat = targ(i, :);
%     normParams_targ(i, :) = [mean(feat), std(feat)];
%     targ(i, :) = (feat - normParams_targ(i, 1))/normParams_targ(i, 2);
%     
% end
%     
% figure; hold on;
% for i=1:size(targ,1)
%     subplot(size(targ,1),1,i)
%     h1 = histogram(targ(i,:),100); hold on;
%     h2 = histogram(targnotnorm(i,:),100);
%     legend('normalized','not normalized');
%     title(['Hist of normalized ',targetCells{1,i}]); xlabel([targetCells{1,i},'[units]']);
% end
% hold off;

%% Performance over Hidden neuron Number
% NumOfRepeats = 10;
% HiddenN = [2,5,10,15,20,25,30,35,38,40,42,45,50,55,60];
NumOfRepeats = 10;
HiddenN = [5,10,15,20];
[ NN_Mean_over_HN_num,NN_stdev_over_HN_num ] = NN_Perf_over_HNnum( NumOfRepeats,HiddenN,sampl,targ,1 );

%% Performance over samples quantity
NumOfRepeats = 7;
HiddenN = 10;
dataPointsNum = [10000,30000,50000,80000,100000,120000,180000,220000,320000,420000,...
    520000,591000];
[ NN_Mean_over_sampl_num,NN_stdev_over_sampl_num ] = NN_Perf_over_sampl_num( NumOfRepeats,HiddenN,dataPointsNum,sampl,targ,1 );


%% visualiesd the results in the W matrix
NN_weights_matrix_plot( net,parametersCells )

%% visualiesd the results in the W matrix 4X4 for each neuron

 % multiply the output weights with the neurons outputs
weightsInput = net.IW;  weightsInput = cell2mat(weightsInput);
weightsOutput = net.LW;  weightsOutput = cell2mat(weightsOutput);
bias = net.b;  bias = cell2mat(bias);
weights = horzcat(diag(weightsOutput)*weightsInput); 

[x,y] = meshgrid(1:4,1:4);   %# Create x and y coordinates for the strings
HiddenNumCells = num2cell(1:4);
textStrings = num2str(weights(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
%
figure;
n = 3; % the offset in which the W_ij starts (if we have c_i or not). '3' = no c_i, '7' = with c_i 
for i=1:size(weights,1)
    h=subplot(2,5,i);
    N_i_weights = [0              ,weights(i,n)    ,weights(i,n+1) ,weights(i,n+2);
                   weights(i,n+3) ,0               ,weights(i,n+4) ,weights(i,n+5);
                   weights(i,n+6) ,weights(i,n+7)  ,0              ,weights(i,n+8);
                   weights(i,n+9) ,weights(i,n+10) ,weights(i,n+11),0             ];
    N_i_weights_T = N_i_weights';
    imagesc(abs(N_i_weights));
    colormap(flipud(gray));
    title(['Neuron ',num2str(i), ' weights']);
    start_index = 2*size(weights,1)+i;
    jumps = size(weights,1);
    N_i_textStrings = textStrings(start_index:jumps:end); % get the weights as text
    N_i_textStrings = vertcat(' ',N_i_textStrings(1:3),...
                        N_i_textStrings(4),' ',N_i_textStrings(5:6),...
                        N_i_textStrings(7:8),' ',N_i_textStrings(9),...
                        N_i_textStrings(10:12),' ');
    hStrings = text(y(:),x(:),N_i_textStrings,'HorizontalAlignment','center');
    midValue1 = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(abs(N_i_weights_T(:)) > abs(midValue1),1,3);%# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    set(gca,'XTick',1:4,...     %# Change the axes tick marks
        'XTickLabel',HiddenNumCells,...  %#   and tick labels
        'YTick',1:4,...
        'YTickLabel',HiddenNumCells,...
        'TickLength',[0 0]);
end

MML = MatsuokaML();
figure; hold on;
for i = 1:size(weights,1);
    h1 = subplot(2,5,i);
    MML.drawNet(weights(i,3:14)) 
end
hold off

%% prepare the data points that we want
dataPointsNum = 50000;
[ samplNewSymmetric,TargetNewSymmetric ] = matsuoka_symmetric( dataPointsNum,sampl,targ );

HiddenN = 10;
net = feedforwardnet(HiddenN);
[net, tr] = train(net, samplNewSymmetric, TargetNewSymmetric);
bestNet = net;
bestValidSoFar = tr.best_vperf;
for i=1:10
    net = feedforwardnet(HiddenN);
    [net, tr] = train(net, samplNewSymmetric, TargetNewSymmetric);
    if ((tr.best_vperf > bestValidSoFar) && (tr.best_vperf>tr.best_perf))
        bestNet = net;
        bestValidSoFar = tr.best_vperf;
    end
end
%% rearrenging the CPG in a unique way
dataPointsNum = 50000;
[ samplNewUniq,TargetNewUniq ] = matsuoka_uniq( dataPointsNum,sampl,targ )
%% Train NN with 1 layer
HiddenN = 10;
net = feedforwardnet(HiddenN);
% net = fitnet(HiddenN);
[net, tr] = train(net, sampl, targ);

%% Train NN with 2 layers
HiddenN = 10;
net = feedforwardnet([HiddenN,HiddenN]);
[net, tr] = train(net, sampl, targ);
%% Train NN with 3 layers
HiddenN = 12;
net = feedforwardnet([HiddenN,HiddenN,HiddenN]);
[net, tr] = train(net, sampl, targ);

%% Calculating the R^2:

inputs = sampl;
outputs = targ;
netOut = net(inputs);
err = outputs-netOut;
errVar = var(err,0,2);
inputVar = var(outputs,0,2);

R_squar = 1-(errVar/inputVar);



