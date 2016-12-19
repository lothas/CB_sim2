

clc; close all; clear all;

%% Load varibles:
load('MatsRandomRes_16_12_2016.mat','nSims','results','periods')

%% define the cell names for weight color matrix
% parametersCells = {'\tau','b','c_1','c_2','c_3','c_4',...
%     'w_{12}','w_{13}','w_{14}',...
%     'w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}',...
%     'w_{41}','w_{42}','w_{43}'};
% parametersCells = {'\tau','b',...
%     'w_{12}','w_{13}','w_{14}',...
%     'w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}',...
%     'w_{41}','w_{42}','w_{43}',...
%     'sum W','\Pi W'};
parametersCells = {'\tau','b',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};

%% preprocess the data (get rid off unwnted periods)
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error

ids = find(ids_period & ids_error);

figure;
histogram(max(horzcat(results(ids).perError2)',[],2),1000);
title('Hist of period error'); xlabel('period error');

% wantedParam = [tau ,b ,c1-c4, W_ij]
wantedParam = [1,2,7:18];
% ids = find(periods>0.05 & periods<3);
MatsuokaParam = vertcat(results(ids).seq)';
MatsuokaParam = MatsuokaParam(wantedParam,:);
prod_W = prod(MatsuokaParam,1);
sum_W = sum(MatsuokaParam,1);
sampl=vertcat(MatsuokaParam);%,sum_W,prod_W);
targPeriod = periods(ids);
targFreq = 1./periods(ids);

targ = targFreq; % what is the target? period or frequency?

targnotnorm = targ;
clear results periods

%%
figure;
histogram(targ,100);
title('Hist of period'); xlabel('period [sec]');

figure;
histogram(sampl(1,:),100);
title('Hist of \tau'); xlabel('\tau [sec]');

%% Normalize samples
normParams = zeros(size(sampl, 1), 2);
for i = 1:size(sampl, 1)
    feat = sampl(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    sampl(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    
end

normTarg = [mean(targ), std(targ)];
targ = (targ - normTarg(1, 1))/normTarg(1, 2);
    
figure;
h1 = histogram(targ,100); hold on;
h2 = histogram(targnotnorm,100);
legend('normalized','not normalized');
title('Hist of normalized period'); xlabel('period [sec]');

%% Performance over Hidden neuron Number
NumOfRepeats = 10;
HiddenN = [2,5,10,15,20,25,30,35,38,40,42,45,50,55,60];
netMseTrain = zeros(NumOfRepeats,length(HiddenN));
netMseValidation = zeros(NumOfRepeats,length(HiddenN));
netMseTest = zeros(NumOfRepeats,length(HiddenN));
% net.trainParam.epochs = 300;
% net.trainParam.goal = 1e-5;
for i=1:length(HiddenN)
    for j=1:NumOfRepeats
        net = feedforwardnet(HiddenN(1,i));
        [net, tr] = train(net, sampl, targ);
        netMseTrain(j,i) = tr.best_perf;
        netMseValidation(j,i) = tr.best_vperf;
        netMseTest(j,i) = tr.best_tperf;
        
        clear net tr
    end

end
%
% save('NNdata.mat','netMseTrain','netMseValidation','netMseTest','HiddenN')
meanMseTrain = mean(netMseTrain,1);
meanMseValidation = mean(netMseValidation,1);
meanMseTest = mean(netMseTest,1);

stdMseTrain = std(netMseTrain,0,1);
stdMseValidation = std(netMseValidation,0,1);
stdMseTest = std(netMseTest,0,1);

figure;
errorbar(HiddenN,meanMseTrain,stdMseTrain); hold on
errorbar(HiddenN,meanMseValidation,stdMseValidation);
errorbar(HiddenN,meanMseTest,stdMseTest);
legend('Train group','validation group','Test group');
title('network performance (MSE) over hidden neurons num (target=frequency)');
xlabel('Hidden Neuron Num');
ylabel('MSE');

%% Performance over samples quantity
NumOfRepeats = 10;
HiddenN = 30;
dataPointsNum = [1000,5000,10000,20000,30000,40000,50000,60000,70000,80000,90000];
netMseTrain = zeros(NumOfRepeats,length(dataPointsNum));
netMseValidation = zeros(NumOfRepeats,length(HiddenN));
netMseTest = zeros(NumOfRepeats,length(HiddenN));
% net.trainParam.epochs = 300;
% net.trainParam.goal = 1e-5;
for i=1:length(dataPointsNum)
    for j=1:NumOfRepeats
        net = feedforwardnet(HiddenN);
        dataInd4Train = randsample(length(targ),dataPointsNum(1,i));
        sampl4train = sampl(:,dataInd4Train);
        targ4train = targ(1,dataInd4Train);
        [net, tr] = train(net, sampl4train, targ4train);
        netMseTrain(j,i) = tr.best_perf;
        netMseValidation(j,i) = tr.best_vperf;
        netMseTest(j,i) = tr.best_tperf;
        
        clear net tr dataInd4Train sampl4train targ4train
    end

end
%
save('NNdata_numPoints.mat','netMseTrain','netMseValidation','netMseTest','HiddenN','dataPointsNum')
meanMseTrain = mean(netMseTrain,1);
meanMseValidation = mean(netMseValidation,1);
meanMseTest = mean(netMseTest,1);

stdMseTrain = std(netMseTrain,0,1);
stdMseValidation = std(netMseValidation,0,1);
stdMseTest = std(netMseTest,0,1);

figure;
errorbar(dataPointsNum,meanMseTrain,stdMseTrain); hold on
errorbar(dataPointsNum,meanMseValidation,stdMseValidation);
errorbar(dataPointsNum,meanMseTest,stdMseTest);
legend('Train group','validation group','Test group');
title('network performance over hidden neurons num MSE BEST');
xlabel('data points num');
ylabel('MSE');

%% visualiesd the results in the W matrix

 % multiply the output weights with the neurons outputs
weightsInput = net.IW;  weightsInput = cell2mat(weightsInput);
weightsOutput = net.LW;  weightsOutput = cell2mat(weightsOutput);
bias = net.b;  bias = cell2mat(bias);
weights = horzcat(diag(weightsOutput)*weightsInput); 

%

[x,y] = meshgrid(1:size(weights,2),1:size(weights,1));   %# Create x and y coordinates for the strings
HiddenNumCells = num2cell((1:size(weights,1)));
textStrings = num2str(weights(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding

figure;
imagesc(abs(weights));
colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)
title(['Weights in a NN with',num2str(HiddenN),'hidden neurons']);
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(abs(weights(:)) > midValue,1,3);%# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gca,'XTick',1:size(weights,2),...     %# Change the axes tick marks
        'XTickLabel',parametersCells,...  %#   and tick labels
        'YTick',1:size(weights,1),...
        'YTickLabel',HiddenNumCells,...
        'TickLength',[0 0]);
    

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



