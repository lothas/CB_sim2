
%% prepare the data points that we want
dataPointsNum = 50000;

dataInd4Train = randsample(length(targ),dataPointsNum);
sampl4train = sampl(:,dataInd4Train);
targ4train = targ(1,dataInd4Train);

% sampl4train=[0;0;1;2;3;4;12;13;14;21;23;24;31;32;34;41;42;43];
% targ4train = 1;

%% how many ways we can arrange the Matsuoka network
mixingVectos = perms([1,2,3,4]);


samplNewSymmetric = sampl4train;
TargetNewSymmetric = targ4train;

for j=1:length(mixingVectos)
    
    mixVector = mixingVectos(j,:);
    newSampl = zeros(size(sampl4train,1),dataPointsNum);
    for i=1:dataPointsNum

        newSampl(1:2,i) = sampl4train(1:2,i);
%         newSampl(3:6,i) = sampl4train(mixVector+2,i); % no C

        mat = [0             ,sampl4train(3,i)   ,sampl4train(4,i) ,sampl4train(5,i);
             sampl4train(6,i),0                  ,sampl4train(7,i) ,sampl4train(8,i);
             sampl4train(9,i),sampl4train(10,i)  ,0                ,sampl4train(11,i);
             sampl4train(12,i),sampl4train(13,i) ,sampl4train(14,i),0               ];

        mat = mat(mixVector,mixVector);

        newSampl(3:14,i) = [mat(1,2);mat(1,3);mat(1,4);...
                            mat(2,1);mat(2,3);mat(2,4);...
                            mat(3,1);mat(3,2);mat(3,4);...
                            mat(4,1);mat(4,2);mat(4,3)];
                        
    end

    samplNewSymmetric = horzcat(samplNewSymmetric,newSampl);
    TargetNewSymmetric = horzcat(TargetNewSymmetric,targ4train);
end

%%%% training the NN:
if 0
    HiddenN = 10;
    net = feedforwardnet(HiddenN);
    [net, tr] = train(net, samplNewSymmetric, TargetNewSymmetric);
end

%% making the unique structure:
% samplNewUniq = 
% TargetNewUniq
    
newSampl = zeros(size(sampl4train,1),dataPointsNum);

for i=1:dataPointsNum

    newSampl(1:2,i) = sampl4train(1:2,i);

    mat = [0             ,sampl4train(7,i)  ,sampl4train(8,i) ,sampl4train(9,i);
        sampl4train(10,i),0                 ,sampl4train(11,i),sampl4train(12,i);
        sampl4train(13,i),sampl4train(14,i) ,0                ,sampl4train(15,i);
        sampl4train(16,i),sampl4train(17,i) ,sampl4train(18,i),0               ];
    
    biggestWeights = max(mat,[],1);
    [~,sortedWeightsIND] = sort(biggestWeights,2);
%     newSampl(2+(1:4),i) = sampl4train(sortedWeightsIND+2,i); % no C_i
    
    mat = mat(sortedWeightsIND',sortedWeightsIND');

    newSampl(3:14,i) = [mat(1,2);mat(1,3);mat(1,4);...
                        mat(2,1);mat(2,3);mat(2,4);...
                        mat(3,1);mat(3,2);mat(3,4);...
                        mat(4,1);mat(4,2);mat(4,3)];


end
% 
% samplNewUniq = horzcat(samplNewSymmetric,newSampl);
% TargetNewUniq = horzcat(TargetNewSymmetric,targ4train);
 %%%% training the NN:
if 1
    HiddenN = 10;
    net = feedforwardnet(HiddenN);
    [net, tr] = train(net, newSampl, targ4train);
end   

%% cheking symmetry with uniqeness

    
newSampl = zeros(size(samplNewSymmetric,1),length(TargetNewSymmetric));

for i=1:length(TargetNewSymmetric)

    newSampl(1:2,i) = samplNewSymmetric(1:2,i);

    mat = [0             ,samplNewSymmetric(7,i)  ,samplNewSymmetric(8,i) ,samplNewSymmetric(9,i);
        samplNewSymmetric(10,i),0                 ,samplNewSymmetric(11,i),samplNewSymmetric(12,i);
        samplNewSymmetric(13,i),samplNewSymmetric(14,i) ,0                ,samplNewSymmetric(15,i);
        samplNewSymmetric(16,i),samplNewSymmetric(17,i) ,samplNewSymmetric(18,i),0               ];
    
    biggestWeights = max(mat,[],1);
    [~,sortedWeightsIND] = sort(biggestWeights,2);
    newSampl(3:6,i) = samplNewSymmetric(sortedWeightsIND+2,i);
    
    mat = mat(sortedWeightsIND',sortedWeightsIND');

    newSampl(7:18,i) = [mat(1,2);mat(1,3);mat(1,4);...
                        mat(2,1);mat(2,3);mat(2,4);...
                        mat(3,1);mat(3,2);mat(3,4);...
                        mat(4,1);mat(4,2);mat(4,3)];


end
% 
firstPerm = newSampl(1,1:dataPointsNum);
for i=1:24
    check = sum(newSampl(1,(dataPointsNum*i+1):dataPointsNum*(i+1)) == firstPerm);
    if check ~= dataPointsNum
        disp(['check fail: permutation ' , num2str(i), 'failed']);
    else
        disp(['check ',num2str(i), ' passed'])
    end
end